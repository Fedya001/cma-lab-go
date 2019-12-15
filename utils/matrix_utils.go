package utils

import (
	"cma-lab-go/matrix"
	"math"
	"math/rand"
)

func GenerateSquareMatrix(size int32, min, max float64) *matrix.SquareMatrix {
	data := make([][]float64, 0, size)

	for i := int32(0); i < size; i++ {
		data = append(data, make([]float64, 0, size))
		for j := int32(0); j < size; j++ {
			data[i] = append(data[i], min + rand.Float64() * (max - min))
		}
	}

	// TODO: how to write this in one line ?
	result, _ := matrix.NewSquareMatrix(data)
	return result
}

// Euclid norm is used
func GetNorm(vector []float64) float64 {
	var max float64 = 0
	for _, v := range vector {
		max = math.Max(max, math.Abs(v))
	}
	return max
}

func NormColumn(column *matrix.Column) {
	norm := GetNorm(column.Data)
	for i, _ := range column.Data {
		column.Data[i] /= norm
	}
}

func NormRow(row *matrix.Row) {
	norm := GetNorm(row.Data)
	for i, _ := range row.Data {
		row.Data[i] /= norm
	}
}