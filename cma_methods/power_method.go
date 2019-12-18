package cma_methods

import (
	"cma-lab-go/matrix"
	"cma-lab-go/utils"
	"math"
	"math/cmplx"
	"sync"
)

type EigenvalueCase int

const (
	REAL_EIGENVALUE_CASE EigenvalueCase = iota
	OPPOSITE_PAIRED_EIGENVALUES_CASE
	COMPLEX_EIGENVALUES_CASE
	STUCK_CASE
)

const STOP_THRESHOLD = 1e-11

type Eigenvector struct {
	Value  complex128
	Vector []complex128
}

// Can't use MultiplyMatrixOnColumn func here, because SquareMatrix is real
// and there are no generics in go :(
func ensureEigenvectors(squareMatrix *matrix.SquareMatrix, eigenvectors ...*Eigenvector) bool {
	dim := len(squareMatrix.Data)

	var wg = sync.WaitGroup{}

	for index := 0; index < len(eigenvectors); index++ {
		// simple test on nonzero vector
		allZero := true
		for i := 0; i < dim; i++ {
			if math.Abs(real(eigenvectors[index].Vector[i])) > STOP_THRESHOLD {
				allZero = false
				break
			}
		}

		if allZero {
			return false
		}


		wg.Add(dim)
		var sum float64 = 0
		for i := 0; i < dim; i++ {
			go func(i int) {
				var localSum complex128 = 0
				for j := 0; j < dim; j++ {
					localSum += complex(squareMatrix.Data[i][j], 0) * eigenvectors[index].Vector[j]
				}
				localSum -= eigenvectors[index].Value * eigenvectors[index].Vector[i]
				abs := cmplx.Abs(localSum)
				sum += abs * abs

				wg.Done()
			}(i)
		}
		wg.Wait()

		if math.Sqrt(sum) > STOP_THRESHOLD {
			return false
		}
	}

	return true
}

// REAL_EIGENVALUE_CASE
func makeRealEigenvector(prev, cur []float64) *Eigenvector {
	converted := make([]complex128, 0, len(prev))
	for i, _ := range cur {
		converted = append(converted, complex(cur[i], 0))
	}

	var max float64 = 0
	for i := 0; i < len(cur); i++ {
		if prev[i] != 0 {
			max = math.Max(max, cur[i]/prev[i])
		}
	}

	return &Eigenvector{complex(max, 0), converted}
}

// OPPOSITE_PAIRED_EIGENVALUES_CASE
func makeOppositeEigenvectors(last, prev, cur []float64) []*Eigenvector {
	var max float64 = 0
	for i := 0; i < len(cur); i++ {
		if last[i] != 0 {
			max = math.Max(max, cur[i]/last[i])
		}
	}
	lambda := math.Sqrt(max)

	plusVector := make([]complex128, 0, len(cur))
	minusVector := make([]complex128, 0, len(cur))

	for i, _ := range cur {
		plusVector = append(plusVector, complex(cur[i]+lambda*prev[i], 0))
		minusVector = append(minusVector, complex(cur[i]-lambda*prev[i], 0))
	}

	return []*Eigenvector{
		{Value: complex(lambda, 0), Vector: plusVector},
		{Value: complex(-lambda, 0), Vector: minusVector},
	}
}

// COMPLEX_EIGENVALUES_CASE
func makeComplexEigenvector(start, last, prev, cur []float64) ([]*Eigenvector, bool) {
	dim := len(cur)

	var r float64 = 0
	for i := 0; i < dim; i++ {
		r = math.Max(r, math.Sqrt(math.Abs((last[i]*cur[i]-prev[i]*prev[i])/
			(start[i]*prev[i]-last[i]*last[i]))))
	}

	var cosTetta float64 = 0
	for i := 0; i < dim; i++ {
		c := (cur[i]+r*r*last[i]) / (2*r*prev[i])
		if math.Abs(c) > math.Abs(cosTetta) {
			cosTetta = c
		}
	}

	if math.Abs(cosTetta) > 1 {
		return make([]*Eigenvector, 0), false
	}

	sinTetta := math.Sqrt(1 - cosTetta*cosTetta)

	lambdaOne := complex(r*cosTetta, r*sinTetta)
	lambdaTwo := complex(r*cosTetta, -r*sinTetta)

	// Create eigenvectors
	vectorOne := make([]complex128, 0, dim)
	vectorTwo := make([]complex128, 0, dim)

	for i := 0; i < dim; i++ {
		vectorOne = append(vectorOne, complex(cur[i], 0)-lambdaTwo*complex(prev[i], 0))
		vectorTwo = append(vectorTwo, complex(cur[i], 0)-lambdaOne*complex(prev[i], 0))
	}

	return []*Eigenvector{
		{Value: lambdaOne, Vector: vectorOne},
		{Value: lambdaTwo, Vector: vectorTwo},
	}, true
}

func FindMaxEigenvalues(squareMatrix *matrix.SquareMatrix, initApprox *matrix.Column) ([]*Eigenvector, EigenvalueCase) {
	// middle = iteration before prev, last = iteration before middle

	start, _ := matrix.MultiplyMatrixOnColumn(squareMatrix, initApprox)
	var last, prev, cur *matrix.Column

	// new modification : process a bunch of 4 vectors
	// to simplify evaluations
	for {
		// process 4 iterations at once
		last, _ = matrix.MultiplyMatrixOnColumn(squareMatrix, start)
		prev, _ = matrix.MultiplyMatrixOnColumn(squareMatrix, last)
		cur, _ = matrix.MultiplyMatrixOnColumn(squareMatrix, prev)

		// Evaluate all eigenvalues using 3 formulas (for each case)

		// case 1 : real eigenvalue
		realEigenvector := makeRealEigenvector(prev.Data, cur.Data)
		if ensureEigenvectors(squareMatrix, realEigenvector) {
			return []*Eigenvector{realEigenvector}, REAL_EIGENVALUE_CASE
		}

		// case 2 : opposite pairing eigenvalues
		oppositeEigenvectors := makeOppositeEigenvectors(last.Data, prev.Data, cur.Data)
		if ensureEigenvectors(squareMatrix, oppositeEigenvectors...) {
			return oppositeEigenvectors, OPPOSITE_PAIRED_EIGENVALUES_CASE
		}

		// case 3 : complex eigenvalues case
		complexEigenvectors, ok := makeComplexEigenvector(start.Data, last.Data, prev.Data, cur.Data)
		if ok && ensureEigenvectors(squareMatrix, complexEigenvectors...) {
			return complexEigenvectors, COMPLEX_EIGENVALUES_CASE
		}

		utils.NormColumn(cur)
		start, _ = matrix.MultiplyMatrixOnColumn(squareMatrix, cur)
		utils.NormColumn(start)
	}

}
