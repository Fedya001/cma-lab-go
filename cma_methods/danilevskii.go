package cma_methods

import (
	"cma-lab-go/matrix"
	"math"
	"sync"
)

func multiplyPolynomials(lhs, rhs []float64) []float64 {
	n, m := len(lhs), len(rhs)
	result := make([]float64, n + m - 1)

	for k1, v1 := range lhs {
		for k2, v2 := range rhs {
			result[k1+k2] += v1 * v2
		}
	}

	last := n + m - 2
	for ; last > 0 && math.Abs(result[last]) < 1e-9; {
		last--
	}

	// TODO : check it
	return result[:last + 1]
}

func FindPolynomial(squareMatrixOriginal *matrix.SquareMatrix) []float64 {
	size := len(squareMatrixOriginal.Data)

	// make a deep copy
	squareMatrix := matrix.SquareMatrix{}
	squareMatrix.Data = make([][]float64, size)
	for i := 0; i < size; i++ {
		squareMatrix.Data[i] = make([]float64, size)
		copy(squareMatrix.Data[i], squareMatrixOriginal.Data[i])
	}

	var polynomials [][]float64

	var wg = sync.WaitGroup{}

	prevSlice := size
	for column := size - 2; column >= 0; column-- {

		// swap rows and cols is required
		var maxInd int = column
		for i := 0; i < column; i++ {
			if math.Abs(squareMatrix.Data[column+1][i]) >
				math.Abs(squareMatrix.Data[column+1][maxInd]) {
				maxInd = i
			}
		}

		if maxInd != column {
			for i := 0; i <= column+1; i++ {
				squareMatrix.Data[i][maxInd], squareMatrix.Data[i][column] =
					squareMatrix.Data[i][column], squareMatrix.Data[i][maxInd]
			}

			squareMatrix.Data[maxInd], squareMatrix.Data[column] =
				squareMatrix.Data[column], squareMatrix.Data[maxInd]
		}

		if math.Abs(squareMatrix.Data[column+1][column]) < 1e-9 {
			// split case
			polynom := make([]float64, 0, prevSlice-column-1)
			pol_len := prevSlice - column - 1
			for i := prevSlice - 1; i >= column+1; i-- {
				if pol_len%2 == 0 {
					polynom = append(polynom, -squareMatrix.Data[column+1][i])
				} else {
					polynom = append(polynom, squareMatrix.Data[column+1][i])
				}
			}
			if pol_len%2 == 0 {
				polynom = append(polynom, 1)
			} else {
				polynom = append(polynom, -1)
			}

			prevSlice = column + 1
			polynomials = append(polynomials, polynom)
			continue
		}

		// save the [column + 1] row
		var baseRowBuf matrix.Row
		baseRowBuf.Data = make([]float64, size)
		copy(baseRowBuf.Data, squareMatrix.Data[column+1])

		// M_n
		for i := 0; i <= column+1; i++ {
			squareMatrix.Data[i][column] /= squareMatrix.Data[column+1][column]
		}

		wg.Add(size - 1)
		for j := 0; j < size; j++ {
			if j == column {
				continue
			}

			go func(j, column int) {
				for i := 0; i <= column+1; i++ {
					squareMatrix.Data[i][j] -= squareMatrix.Data[i][column] * squareMatrix.Data[column+1][j]
				}

				wg.Done()
			}(j, column)
		}
		wg.Wait()

		// M_n^{-1}
		buf, _ := matrix.MultiplyRowOnMatrix(&baseRowBuf, &squareMatrix)

		copy(squareMatrix.Data[column], buf.Data)
	}

	// last polynomial
	polynom := make([]float64, 0, prevSlice)
	for i := prevSlice - 1; i >= 0; i-- {
		if prevSlice%2 == 0 {
			polynom = append(polynom, -squareMatrix.Data[0][i])
		} else {
			polynom = append(polynom, squareMatrix.Data[0][i])
		}
	}

	if prevSlice%2 == 0 {
		polynom = append(polynom, 1)
	} else {
		polynom = append(polynom, -1)

	}

	polynomials = append(polynomials, polynom)

	// multiply polynomials
	result := make([]float64, 0)
	result = append(result, 1)
	for i, _ := range polynomials {
		result = multiplyPolynomials(result, polynomials[i])
	}

	return result
}
