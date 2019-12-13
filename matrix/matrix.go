package matrix

import (
	"fmt"
	"sync"
)

// TODO: I don't know whether it should be global
var wg = sync.WaitGroup{}

type Column struct {
	data []float64
}

type Row struct {
	data []float64
}

type SquareMatrix struct {
	data [][]float64
}

// All slices are copied by pointer, therefore
// be careful with side effects

// Initialization functions
func NewColumn(data []float64) *Column {
	return &Column{data: data}
}

func NewRow(data []float64) *Row {
	return &Row{data: data}
}

func NewSquareMatrix(data [][]float64) (*SquareMatrix, error) {
	for i, _ := range data {
		if len(data[i]) != len(data) {
			return &SquareMatrix{}, fmt.Errorf("inconsistent matrix size")
		}
	}
	return &SquareMatrix{data: data}, nil
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
	if len(matrix.data) != len(column.data) {
		return &Column{}, fmt.Errorf("inconsistent matrix and column sizes")
	}

	var result Column
	result.data = make([]float64, len(matrix.data))

	wg.Add(len(result.data))
	for i := 0; i < len(result.data); i++ {
		go func(i int) {
			result.data[i] = multiplyVectors(matrix.data[i], column.data)
			wg.Done()
		}(i)
	}
	wg.Wait()

	return &result, nil
}

func MultiplyRowOnMatrix(row *Row, matrix *SquareMatrix) (*Row, error) {
	if len(row.data) != len(matrix.data) {
		return &Row{}, fmt.Errorf("inconsistent row and matrix sizes")
	}

	var result Row
	result.data = make([]float64, len(matrix.data))

	wg.Add(len(result.data))
	for i := 0; i < len(result.data); i++ {
		go func(i int) {
			result.data[i] = multiplyVectors(row.data, matrix.data[i])
			wg.Done()
		}(i)
	}
	wg.Wait()

	return &result, nil
}

// TODO: how to remove all this copy-paste blocks ?
func MultiplyMatrices(lhs, rhs *SquareMatrix) (*SquareMatrix, error) {
	if len(lhs.data) != len(rhs.data) {
		return &SquareMatrix{}, fmt.Errorf("inconsistent matrices sizes")
	}

	var result SquareMatrix
	result.data = make([][]float64, len(lhs.data))

	wg.Add(len(lhs.data) * len(rhs.data)) // a green thread per each pair of vectors
	for i := 0; i < len(lhs.data); i++ {
		result.data[i] = make([]float64, len(lhs.data))
		for j := 0; j < len(rhs.data); j++ {
			go func(i, j int) {
				result.data[i][j] = multiplyVectors(lhs.data[i], rhs.data[j])
				wg.Done()
			}(i, j)
		}
	}
	wg.Wait()

	return &result, nil
}
