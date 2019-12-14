package io

import (
	"bufio"
	"cma-lab-go/matrix"
	"io"
	"os"
	"strconv"
)

type FileMatrixReader struct {
	file    *os.File
	scanner *bufio.Scanner
}

func NewFileMatrixReader(filename string) (*FileMatrixReader, error) {
	matrixReader := FileMatrixReader{}
	file, err := os.Open(filename)
	if err != nil {
		return &matrixReader, err
	}

	matrixReader.file = file
	matrixReader.scanner = bufio.NewScanner(file)
	matrixReader.scanner.Split(bufio.ScanWords)
	return &matrixReader, nil
}

func (reader *FileMatrixReader) readVector(size int) ([]float64, error) {
	result := make([]float64, 0, size)
	for i := 0; i < size; i++ {
		reader.scanner.Scan()
		x, err := strconv.ParseFloat(reader.scanner.Text(), 64)
		if err != nil {
			return result, err
		}
		result = append(result, x)
	}
	return result, nil
}

func (reader *FileMatrixReader) ReadDimension() (int, error) {
	reader.scanner.Scan()
	return strconv.Atoi(reader.scanner.Text())
}

func (reader *FileMatrixReader) ReadRow(size int) (*matrix.Row, error) {
	row, err := reader.readVector(size)
	return matrix.NewRow(row), err
}

func (reader *FileMatrixReader) ReadColumn(size int) (*matrix.Column, error) {
	col, err := reader.readVector(size)
	return matrix.NewColumn(col), err
}

func (reader *FileMatrixReader) ReadMatrix(size int) (*matrix.SquareMatrix, error) {
	result := make([][]float64, 0, size)
	for i := 0; i < size; i++ {
		row, err := reader.readVector(size)
		if err != nil {
			return &matrix.SquareMatrix{}, err
		}
		result = append(result, row)
	}

	return matrix.NewSquareMatrix(result)
}

func (reader *FileMatrixReader) Close() error {
	return reader.file.Close()
}

type MatrixReader interface {
	io.Closer
	ReadDimension() (int, error)
	ReadMatrix(size int) (*matrix.SquareMatrix, error)
	ReadRow(size int) (*matrix.Row, error)
	ReadColumn(size int) (*matrix.Column, error)
}
