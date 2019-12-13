package utils

import (
	"cma-lab-go/matrix"
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