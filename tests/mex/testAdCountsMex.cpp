#include "mex.h"
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "kalmanLogLikelihoodAd.h"

namespace {
void require(bool condition, const std::string& name, const char* message)
{
    if (!condition)
        throw std::runtime_error(name + ": " + message);
}

void runCase(const std::string& name, std::size_t m, std::size_t n,
             const std::vector<std::size_t>& offsets,
             const std::vector<std::size_t>& columns,
             const std::vector<double>& values)
{
    SparseMatrix<double> A(m, n);
    A.row().resize(m + 1);
    A.col().resize(values.size());
    A.val().resize(values.size());
    for (std::size_t i = 0; i <= m; ++i)
        A.row()(i) = offsets[i];
    for (std::size_t k = 0; k < values.size(); ++k) {
        A.col()(k) = columns[k];
        A.val()(k) = values[k];
    }

    vector<double> weights(n);
    vector<std::size_t> weightIndices(n);
    matrix<double> output(m, n);
    matrix<std::size_t> outputIndices(m, n);
    for (std::size_t j = 0; j < n; ++j) {
        weights(j) = static_cast<double>(j + 1);
        weightIndices(j) = j;
    }
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j)
            outputIndices(i, j) = n + i * n + j;

    AD_As_bdx<double> operation(&A, &weights, &output,
                               &weightIndices, &outputIndices);
    std::size_t countedRows = 0, countedData = 0;
    operation.count(countedRows, countedData);
    require(countedRows == m * n, name, "incorrect row count");
    require(countedData == values.size(), name, "incorrect derivative count");

    const std::size_t guard = std::numeric_limits<std::size_t>::max();
    vector<double> derivatives(countedData + 1, -12345.0);
    vector<std::size_t> indices(countedData + 1, guard);
    vector<std::size_t> rows(countedRows + 1, guard);
    vector<std::size_t> resultIndices(countedRows + 1, guard);
    std::size_t nr = 0, nd = 0;
    operation.fAD(derivatives, indices, rows, resultIndices, nr, nd);
    require(nr == countedRows && nd == countedData, name,
            "recorded tape does not match count()");
    require(derivatives(countedData) == -12345.0 &&
            indices(countedData) == guard &&
            resultIndices(countedRows) == guard && rows(countedRows) == guard,
            name, "write past the counted capacity");
    rows(nr) = nd;

    std::size_t entry = 0;
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t row = i * n + j;
            require(resultIndices(row) == outputIndices(i, j) &&
                    rows(row) == entry, name, "incorrect tape row mapping");
            double expected = 0.0;
            for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
                if (columns[k] == j) {
                    expected += values[k] * weights(j);
                    require(derivatives(entry) == values[k] && indices(entry) == j,
                            name, "incorrect derivative value or column");
                    ++entry;
                }
            }
            require(output(i, j) == expected && rows(row + 1) == entry,
                    name, "incorrect output or tape row length");
        }
    }
    require(entry == countedData, name, "not all stored entries were tested");
}
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray*[])
{
    if (nrhs != 0 || nlhs > 1)
        mexErrMsgIdAndTxt("kalmanrates:AdCountTestArguments", "Use testAdCountsMex() without inputs.");
    try {
        runCase("square", 3, 3, {0, 2, 3, 4}, {0, 2, 1, 0}, {1, -2, 3, 4});
        runCase("rectangular", 2, 4, {0, 2, 3}, {1, 3, 0}, {2, 4, -3});
        runCase("empty row", 3, 2, {0, 1, 1, 2}, {0, 1}, {1, 5});
        runCase("stored zero", 2, 3, {0, 2, 3}, {0, 2, 1}, {0, 7, -4});
        runCase("empty matrix", 2, 3, {0, 0, 0}, {}, {});
        runCase("single column", 2, 1, {0, 1, 2}, {0, 0}, {2, -1});
        if (nlhs == 1)
            plhs[0] = mxCreateDoubleScalar(6);
    } catch (const ExcKal& error) {
        mexErrMsgIdAndTxt("kalmanrates:AdCountTestFailure", "%s", error.msg_.c_str());
    } catch (const std::exception& error) {
        mexErrMsgIdAndTxt("kalmanrates:AdCountTestFailure", "%s", error.what());
    }
}
