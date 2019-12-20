package io

import (
	"cma-lab-go/matrix"
	"fmt"
)

type ConsoleMatrixWriter struct {
}

func NewConsoleMatrixWriter() *ConsoleMatrixWriter {
	return &ConsoleMatrixWriter{}
}

func (writer *ConsoleMatrixWriter) writeSlice(delimiter string, values []float64) {
	for _, v := range values {
		fmt.Printf("\t%v" + delimiter, v)
	}
}

func (writer *ConsoleMatrixWriter) WriteRow(row *matrix.Row) {
	writer.writeSlice(" ", row.Data)
}

func (writer *ConsoleMatrixWriter) WriteColumn(column *matrix.Column) {
	writer.writeSlice("\n", column.Data)
}

func (writer *ConsoleMatrixWriter) WriteMatrix(matrix *matrix.SquareMatrix) {
	for i := 0; i < len(matrix.Data); i++ {
		writer.writeSlice(" ", matrix.Data[i])
		fmt.Printf("\n")
	}
}

type MatrixWriter interface {
	WriteMatrix(*matrix.SquareMatrix)
	WriteRow(row *matrix.Row)
	WriteColumn(column *matrix.Column)
}
