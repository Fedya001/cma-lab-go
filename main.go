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
	matrices := []*matrix.SquareMatrix{}
	for _, filename := range []string{"data/sampleA.txt", "data/sampleB.txt", "data/matrixA.txt", "data/matrixB.txt", } {
		matrixReader, _ := io.NewFileMatrixReader(filename)
		dim, _ := matrixReader.ReadDimension()
		m, _ := matrixReader.ReadMatrix(dim)
		matrices = append(matrices, m)
	}

	var matrixWriter io.MatrixWriter = io.NewConsoleMatrixWriter()

	// 1. Power method
	fmt.Println("1. Power method")
	for _, m := range matrices {
		matrixWriter.WriteMatrix(m)

		initApprox := make([]float64, len(m.Data))
		initApprox[0] = 1

		eigenvectors, methodCase, iterations := cma_methods.FindMaxEigenvalues(m, matrix.NewColumn(initApprox))
		fmt.Printf("Iterations count = %v\n\n", iterations)

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

	// 2. Danilevskii method
	fmt.Println("\n2. Danilevskii method")
	for _, m := range matrices {
		matrixWriter.WriteMatrix(m)

		polynomial := cma_methods.FindPolynomial(m)

		first := true
		for i, v := range polynomial {
			if !first {
				fmt.Print(" + ")
			}
			fmt.Printf("%.8fx^{%v}", v, i)
			first = false
		}
		fmt.Printf("\n\n")
		fmt.Printf("Eigenvalues = %v\n\n", cma_methods.FindPolynomialRoots(polynomial))
	}

	// 3. QR-algorithm
	fmt.Println("3. QR-algorithm")
	for _, m := range matrices {
		matrixWriter.WriteMatrix(m)

		eigenvalues, eigenvectors, iterations := cma_methods.SolveQR(m)
		fmt.Printf("Iterations count = %v\n\n", iterations)
		fmt.Println(eigenvalues)

		for i, _ := range eigenvectors {
			fmt.Printf("Eigenvalue = %v\n", eigenvalues[i])
			fmt.Printf("Eigenvector = %v\n", eigenvectors[i])
		}
		fmt.Println()
	}
}
