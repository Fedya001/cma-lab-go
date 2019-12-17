package main

import (
	"cma-lab-go/cma_methods"
	"cma-lab-go/io"
	"cma-lab-go/matrix"
	"fmt"
	"math/rand"
	"runtime"
	"time"
)

func main() {
	fmt.Println("Cma lab on go! Here we go!")
	runtime.GOMAXPROCS(20) // number of CPUs [configure it later]
	rand.Seed(time.Now().UnixNano())

	// Read matrices
	matrixReaderA, _ := io.NewFileMatrixReader("data/matrixA.txt")
	dimA, _ := matrixReaderA.ReadDimension()
	matrixA, _ := matrixReaderA.ReadMatrix(dimA)

	matrixReaderB, _ := io.NewFileMatrixReader("data/matrixB.txt")
	dimB, _ := matrixReaderB.ReadDimension()
	matrixB, _ := matrixReaderB.ReadMatrix(dimB)

	matrices := []*matrix.SquareMatrix{matrixA, matrixB}

	// 1. Power method
	fmt.Println("1. Power method")
	for _, m := range matrices {
		var matrixWriter io.MatrixWriter = io.NewConsoleMatrixWriter()
		matrixWriter.WriteMatrix(m)

		initApprox := make([]float64, len(m.Data))
		initApprox[0] = 1

		eigenvectors, methodCase := cma_methods.FindMaxEigenvalues(m, matrix.NewColumn(initApprox))

		switch methodCase {
		case cma_methods.REAL_EIGENVALUE_CASE:
			fmt.Println("REAL_EIGENVALUE_CASE")
		case cma_methods.OPPOSITE_PAIRED_EIGENVALUES_CASE:
			fmt.Println("OPPOSITE_PAIRED_EIGENVALUES_CASE")
		case cma_methods.COMPLEX_EIGENVALUES_CASE:
			fmt.Println("COMPLEX_EIGENVALUES_CASE")
		case cma_methods.STUCK_CASE:
			fmt.Println("Method got stuck")
			continue
		}

		for _, v := range eigenvectors {
			fmt.Printf("Eigenvalue = %v\n", v.Value)
			fmt.Printf("Eigenvector = %v\n", v.Vector)
		}
	}
}
