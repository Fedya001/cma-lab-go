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

	var matrixReader io.MatrixReader
	matrices := []string{"data/matrixA.txt", "data/matrixB.txt"}

	// power method
	for _, path := range matrices {
		matrixReader, _ = io.NewFileMatrixReader(path)
		dim, _ := matrixReader.ReadDimension()
		m, _ := matrixReader.ReadMatrix(dim)

		var matrixWriter io.MatrixWriter = io.NewConsoleMatrixWriter()
		matrixWriter.WriteMatrix(m)

		initApprox := make([]float64, dim)
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
