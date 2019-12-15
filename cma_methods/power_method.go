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

const STOP_THRESHOLD = 1e-10

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
func makeRealEigenvector(curV, prevU []float64) *Eigenvector {
	converted := make([]complex128, 0, len(prevU))
	for i, _ := range prevU {
		converted = append(converted, complex(prevU[i], 0))
	}

	var max float64 = 0
	for i := 0; i < len(curV); i++ {
		if prevU[i] != 0 {
			max = math.Max(max, curV[i] / prevU[i])
		}
	}

	return &Eigenvector{complex(max, 0), converted}
}

// OPPOSITE_PAIRED_EIGENVALUES_CASE
func makeOppositeEigenvectors(curV, prevV, lastU []float64) []*Eigenvector {
	var max float64 = 0
	for i := 0; i < len(curV); i++ {
		if lastU[i] != 0 {
			max = math.Max(max, curV[i] * utils.GetNorm(prevV) / lastU[i])
		}
	}
	lambda := math.Sqrt(max)

	plusVector := make([]complex128, 0, len(prevV))
	minusVector := make([]complex128, 0, len(prevV))

	for i, _ := range prevV {
		plusVector = append(plusVector, complex(prevV[i]+lambda*lastU[i], 0))
		minusVector = append(minusVector, complex(prevV[i]-lambda*lastU[i], 0))
	}

	return []*Eigenvector{
		{Value: complex(lambda, 0), Vector: plusVector},
		{Value: complex(-lambda, 0), Vector: minusVector},
	}
}

// COMPLEX_EIGENVALUES_CASE
func makeComplexEigenvector(lastU, lastV, middleU, middleV, prevU, prevV, curU, curV []float64) ([]*Eigenvector, bool) {
	dim := len(curV)

	middleVNorm := utils.GetNorm(middleV)
	prevVNorm := utils.GetNorm(prevV)

	var maxInd = 0
	var maxR float64 = 0
	for i := 0; i < dim; i++ {
		denominator := lastU[i]*prevV[i] - middleU[i]*middleU[i]*middleVNorm
		if math.Abs(denominator) > STOP_THRESHOLD { // just != 0
			candidate := (middleV[i]*curV[i]*prevVNorm-prevV[i]*prevV[i]*middleVNorm)/denominator
			if candidate > maxR {
				maxR = candidate
				maxInd = i
			}
		}
	}
	r := math.Sqrt(maxR)

	cosTetta := (curV[maxInd]*prevVNorm + maxR*middleU[maxInd]) /
		(2 * r * prevV[maxInd])

	if math.Abs(cosTetta) > 1 {
		return make([]*Eigenvector, 0), false
	}

	sinTetta := math.Sqrt(1 - cosTetta*cosTetta)

	lambdaOne := complex(r*cosTetta, r*sinTetta)
	lambdaTwo := complex(r*cosTetta, -r*sinTetta)

	// Create eigenvectors
	vectorOne := make([]complex128, 0, len(prevV))
	vectorTwo := make([]complex128, 0, len(prevV))

	for i := 0; i < dim; i++ {
		vectorOne = append(vectorOne, complex(prevV[i], 0)-lambdaTwo*complex(middleU[i], 0))
		vectorTwo = append(vectorTwo, complex(prevV[i], 0)-lambdaOne*complex(middleU[i], 0))
	}

	return []*Eigenvector{
		{Value: lambdaOne, Vector: vectorOne},
		{Value: lambdaTwo, Vector: vectorTwo},
	}, true
}

func doIteration(squareMatrix *matrix.SquareMatrix, curU, curV *matrix.Column) (*matrix.Column, *matrix.Column) {
	nextU := matrix.NewColumn(make([]float64, len(squareMatrix.Data)))
	copy(nextU.Data, curV.Data)

	utils.NormColumn(nextU)
	nextV, _ := matrix.MultiplyMatrixOnColumn(squareMatrix, nextU)
	return nextU, nextV
}

func FindMaxEigenvalues(squareMatrix *matrix.SquareMatrix, initApprox *matrix.Column) ([]*Eigenvector, EigenvalueCase) {
	// middle = iteration before prev, last = iteration before middle

	lastU := initApprox
	lastV, _ := matrix.MultiplyMatrixOnColumn(squareMatrix, lastU)

	middleU, middleV := doIteration(squareMatrix, lastU, lastV)
	prevU, prevV := doIteration(squareMatrix, middleU, middleV)

	iterationLimit := 100000

	for {
		curU, curV := doIteration(squareMatrix, prevU, prevV)

		// Evaluate all eigenvalues using 3 formulas (for each case)

		// case 1 : real eigenvalue
		realEigenvector := makeRealEigenvector(curV.Data, prevU.Data)
		if ensureEigenvectors(squareMatrix, realEigenvector) {
			return []*Eigenvector{realEigenvector}, REAL_EIGENVALUE_CASE
		}

		// case 2 : opposite pairing eigenvalues
		oppositeEigenvectors := makeOppositeEigenvectors(curV.Data, prevV.Data, lastU.Data)
		if ensureEigenvectors(squareMatrix, oppositeEigenvectors...) {
			return oppositeEigenvectors, OPPOSITE_PAIRED_EIGENVALUES_CASE
		}

		// case 3 : complex eigenvalues case
		complexEigenvectors, ok := makeComplexEigenvector(
			lastU.Data, lastV.Data, middleU.Data, middleV.Data,
			prevU.Data, prevV.Data, curU.Data, curV.Data)

		if ok && ensureEigenvectors(squareMatrix, complexEigenvectors...) {
			return complexEigenvectors, COMPLEX_EIGENVALUES_CASE
		}

		lastU, lastV = middleU, middleV
		middleU, middleV = prevU, prevV
		prevU, prevV = curU, curV

		iterationLimit--
		if iterationLimit == 0 {
			return []*Eigenvector{}, STUCK_CASE
		}
	}

}
