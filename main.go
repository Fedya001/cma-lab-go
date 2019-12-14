package main

import (
	"cma-lab-go/io"
	"fmt"
	"math/rand"
	"runtime"
	"time"
)

func main() {
	fmt.Println("Cma lab on go! Here we go!")
	runtime.GOMAXPROCS(20) // number of CPUs [configure it later]
	rand.Seed(time.Now().UnixNano())

	const MATRIX_A_PATH = "data/matrixA.txt"

	var matrixReader io.MatrixReader
	matrixReader, _ = io.NewFileMatrixReader(MATRIX_A_PATH)

	dim, _ := matrixReader.ReadDimension()
	m, _ := matrixReader.ReadMatrix(dim)

	var matrixWriter io.MatrixWriter = io.NewConsoleMatrixWriter()
	matrixWriter.WriteMatrix(m)
}
