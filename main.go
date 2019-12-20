package main

import (
	"cma-lab-go/cma_methods"
	"cma-lab-go/io"
	"cma-lab-go/matrix"
	"cma-lab-go/utils"
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

		polynomial, transform, blocksSplit := cma_methods.FindPolynomial(m)

		first := true
		for i, v := range polynomial {
			if !first {
				fmt.Print(" + ")
			}
			fmt.Printf("%vx^{%v}", v, i)
			first = false
		}
		fmt.Printf("\n\n")

		eigenvalues := cma_methods.FindPolynomialRoots(polynomial)
		fmt.Printf("Eigenvalues = %v\n", eigenvalues)

		if blocksSplit {
			fmt.Printf("Unable to recover eigenvectors, because there was a block split\n\n")
		} else {
			for _, eigenvalue := range eigenvalues {
				fmt.Println("")
				eigenvector := make([]float64, len(transform.Data))
				var lambda float64 = 1
				for i := len(transform.Data) - 1; i >= 0; i-- {
					eigenvector[i] = lambda
					lambda *= eigenvalue
				}

				eigenvectorTransformed, _ := matrix.MultiplyMatrixOnColumn(transform, matrix.NewColumn(eigenvector))

				fmt.Printf("Eigenvalue = %v\n", eigenvalue)
				fmt.Printf("Eigenvector = %v\n", eigenvectorTransformed.Data)
			}
		}
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
			fmt.Printf("Eigenvector = %v\n\n", eigenvectors[i])
		}
		fmt.Println()
	}

	// measure QR-algorithm time
	min, max := -1000000000.0, 1000000000.0
	for _, size := range []int{
		10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110,
		120, 130, 140, 150, 160, 170, 180, 190, 200,
	} {
		m := utils.GenerateSquareMatrix(size, min, max)

		start := time.Now()
		_, _, _ = cma_methods.SolveQR(m)
		fmt.Printf("Size = %v => time = %v milliseconds\n", size,
			time.Since(start).Milliseconds())
	}
}
