package matrix

import (
	"fmt"
	"sync"
)

// TODO: I don't know whether it should be global
var wg = sync.WaitGroup{}

type Column struct {
	Data []float64
}

type Row struct {
	Data []float64
}

type SquareMatrix struct {
	Data [][]float64
}

// All slices are copied by pointer, therefore
// be careful with side effects

// Initialization functions
func NewColumn(data []float64) *Column {
	return &Column{Data: data}
}

func NewRow(data []float64) *Row {
	return &Row{Data: data}
}

func NewSquareMatrix(data [][]float64) (*SquareMatrix, error) {
	for i, _ := range data {
		if len(data[i]) != len(data) {
			return &SquareMatrix{}, fmt.Errorf("inconsistent matrix size")
		}
	}
	return &SquareMatrix{Data: data}, nil
}

// Multiplication functions
func multiplyVectors(lhs, rhs []float64) float64 {
	var sum float64 = 0
	for i, _ := range lhs {
		sum += lhs[i] * rhs[i]
	}
	return sum
}

func MultiplyMatrixOnColumn(matrix *SquareMatrix, column *Column) (*Column, error) {
	if len(matrix.Data) != len(column.Data) {
		return &Column{}, fmt.Errorf("inconsistent matrix and column sizes")
	}

	var result Column
	result.Data = make([]float64, len(matrix.Data))

	wg.Add(len(result.Data))
	for i := 0; i < len(result.Data); i++ {
		go func(i int) {
			result.Data[i] = multiplyVectors(matrix.Data[i], column.Data)
			wg.Done()
		}(i)
	}
	wg.Wait()

	return &result, nil
}

func MultiplyRowOnMatrix(row *Row, matrix *SquareMatrix) (*Row, error) {
	if len(row.Data) != len(matrix.Data) {
		return &Row{}, fmt.Errorf("inconsistent row and matrix sizes")
	}

	var result Row
	result.Data = make([]float64, len(matrix.Data))

	wg.Add(len(result.Data))
	for i := 0; i < len(result.Data); i++ {
		go func(i int) {
			result.Data[i] = 0
			for j := 0; j < len(result.Data); j++ {
				result.Data[i] += row.Data[j] * matrix.Data[j][i]
			}
			wg.Done()
		}(i)
	}
	wg.Wait()

	return &result, nil
}

// TODO: how to remove all this copy-paste blocks ?
func MultiplyMatrices(lhs, rhs *SquareMatrix) (*SquareMatrix, error) {
	if len(lhs.Data) != len(rhs.Data) {
		return &SquareMatrix{}, fmt.Errorf("inconsistent matrices sizes")
	}

	var result SquareMatrix
	result.Data = make([][]float64, len(lhs.Data))

	wg.Add(len(lhs.Data) * len(rhs.Data)) // a green thread per each pair of vectors
	for i := 0; i < len(lhs.Data); i++ {
		result.Data[i] = make([]float64, len(lhs.Data))
		for j := 0; j < len(rhs.Data); j++ {
			go func(i, j int) {
				for k := 0; k < len(lhs.Data); k++ {
					result.Data[i][j] += lhs.Data[i][k] * rhs.Data[k][j]
				}
				wg.Done()
			}(i, j)
		}
	}
	wg.Wait()

	return &result, nil
}

func (squareMatrix *SquareMatrix) Transpose() {
	size := len(squareMatrix.Data)

	for i := 1; i < size; i++ {
		for j := 0; j < i; j++ {
			squareMatrix.Data[i][j], squareMatrix.Data[j][i] =
				squareMatrix.Data[j][i], squareMatrix.Data[i][j]
		}
	}
}
