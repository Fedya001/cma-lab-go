package cma_methods

import (
	"cma-lab-go/matrix"
	"cma-lab-go/utils"
	"math"
	"math/cmplx"
)

// rotation rules
// i-th row : cos -sin
// j-th row : sin cos
func rotateLeft(m *matrix.SquareMatrix, i, j int, cos, sin float64) {
	if i >= j {
		panic("You are not supposed to call rotateLeft() with i >= j")
	}

	size := len(m.Data)
	for row := 0; row < size; row++ {
		m.Data[i][row], m.Data[j][row] =
			m.Data[i][row]*cos-m.Data[j][row]*sin,
			m.Data[i][row]*sin+m.Data[j][row]*cos
	}
}

func rotateRight(m *matrix.SquareMatrix, i, j int, cos, sin float64) {
	if i >= j {
		panic("You are not supposed to call rotateRight() with i >= j")
	}

	size := len(m.Data)
	for row := 0; row < size; row++ {
		m.Data[row][i], m.Data[row][j] =
			m.Data[row][i]*cos+m.Data[row][j]*sin,
			-m.Data[row][i]*sin+m.Data[row][j]*cos
	}
}

func directZeroElement(transform, squareMatrix *matrix.SquareMatrix, i, j int) {
	if i <= j {
		panic("You are not supposed to call directZeroElement() with i <= j")
	}

	if squareMatrix.Data[i][j] == 0 {
		panic("directZeroElement() is called with squareMatrix.Data[i][j] == 0")
	}

	j++
	denom := math.Sqrt(squareMatrix.Data[i][j-1]*squareMatrix.Data[i][j-1] +
		squareMatrix.Data[j][j-1]*squareMatrix.Data[j][j-1])
	cos, sin := squareMatrix.Data[j][j-1]/denom, squareMatrix.Data[i][j-1]/denom

	// Q^T
	// j : cos -sin
	// i : sin cos
	rotateLeft(squareMatrix, j, i, cos, -sin)
	rotateRight(squareMatrix, j, i, cos, sin)
}

func doQRIteration(transform, squareMatrix *matrix.SquareMatrix) {
	size := len(squareMatrix.Data)

	// this is crutch for degenerate rotations
	// postponedRightRotationsValues / postponedRightRotationsIndices
	var postponedRRV [][]float64
	var postponedRRI [][]int
	for i := 1; i < size; i++ {
		j := i - 1

		if math.Abs(squareMatrix.Data[i][j]) < STOP_THRESHOLD {
			continue
		}

		denom := math.Sqrt(squareMatrix.Data[j][j] * squareMatrix.Data[j][j] +
			squareMatrix.Data[i][j] * squareMatrix.Data[i][j])
		cos, sin := squareMatrix.Data[j][j]/denom, -squareMatrix.Data[i][j]/denom

		rotateLeft(squareMatrix, j, i, cos, sin)

		postponedRRV = append(postponedRRV, []float64{cos, -sin})
		postponedRRI = append(postponedRRI, []int{j, i})
	}

	for i, _ := range postponedRRI {
		rotateRight(squareMatrix, postponedRRI[i][0], postponedRRI[i][1],
			postponedRRV[i][0], postponedRRV[i][1]) // this is for preserving Q * A * Q^T
		rotateRight(transform, postponedRRI[i][0], postponedRRI[i][1],
			postponedRRV[i][0], postponedRRV[i][1])
	}
}

func extractEigenvalues(squareMatrix *matrix.SquareMatrix) ([]complex128, bool) {
	size := len(squareMatrix.Data)

	var eigenvalues []complex128
	for i := 1; i < size; i++ {
		j := i - 1
		if math.Abs(squareMatrix.Data[i][j]) > STOP_THRESHOLD {
			// block
			trace := squareMatrix.Data[j][j] + squareMatrix.Data[j+1][j+1]
			det := squareMatrix.Data[j+1][j+1]*squareMatrix.Data[j][j] -
				squareMatrix.Data[i][j]*squareMatrix.Data[j][i]

			// l^2 - trace + det = 0
			d := complex(trace*trace-4*det, 0)

			// this is real
			if trace*trace-4*det >= 0 {
				return make([]complex128, 0), false
			}

			eigenvalues = append(eigenvalues,
				(complex(trace, 0)+cmplx.Sqrt(d))/2,
				(complex(trace, 0)-cmplx.Sqrt(d))/2)
			i++
		} else {
			eigenvalues = append(eigenvalues, complex(squareMatrix.Data[j][j], 0))
		}
	}

	// extra check for the real matrix[size - 1][size - 1] element
	// [complex case is processed in cycle]

	if math.Abs(squareMatrix.Data[size-1][size-2]) < STOP_THRESHOLD {
		eigenvalues = append(eigenvalues, complex(squareMatrix.Data[size-1][size-1], 0))
	}

	return eigenvalues, true
}

// fills with zero complex numbers
func extractEigenvectors(squareMatrix, transform *matrix.SquareMatrix) [][]float64 {
	size := len(squareMatrix.Data)

	eigenvectors := make([][]float64, 0, size)
	for i := 1; i < size; i++ {
		j := i - 1
		if math.Abs(squareMatrix.Data[i][j]) > STOP_THRESHOLD {
			eigenvectors = append(eigenvectors, nil, nil)
			i++
		} else {
			vector := make([]float64, 0, size)
			for k := 0; k < size; k++ {
				vector = append(vector, transform.Data[k][j])
			}
			eigenvectors = append(eigenvectors, vector)
		}
	}

	if math.Abs(squareMatrix.Data[size-1][size-2]) < STOP_THRESHOLD {
		vector := make([]float64, 0, size)
		for k := 0; k < size; k++ {
			vector = append(vector, transform.Data[k][size-1])
		}
		eigenvectors = append(eigenvectors, vector)
	}
	return eigenvectors
}

// returns a slice of eigenvalues and slice of eigenvector
// if eigenvalues is complex (Re != 0) => two vectors at this index == nil
func SolveQR(squareMatrixOriginal *matrix.SquareMatrix) ([]complex128, [][]float64, int) {
	size := len(squareMatrixOriginal.Data)

	// deep matrix copy
	squareMatrix := matrix.SquareMatrix{}
	squareMatrix.Data = make([][]float64, size)
	for i := 0; i < size; i++ {
		squareMatrix.Data[i] = make([]float64, size)
		copy(squareMatrix.Data[i], squareMatrixOriginal.Data[i])
	}

	// hessenberg
	// direct rotations
	transform := utils.MakeIdentity(size)
	for j := 0; j < size - 2; j++ {
		for i := j + 2; i < size; i++ {
			if squareMatrix.Data[i][j] != 0 {
				directZeroElement(transform, &squareMatrix, i, j)
			}
		}
	}

	// iterations of QR using rotations (non-direct)
	var prevEigenvalues []complex128
	prevSuccess := false

	iterations := 0
	for {
		doQRIteration(transform, &squareMatrix)
		iterations++

		prev := -2
		curSuccess := true
		for i := 1; i < size; i++ {
			j := i - 1
			if math.Abs(squareMatrix.Data[i][j]) > STOP_THRESHOLD {
				if j == prev+1 {
					// two consequent blocks case
					curSuccess = false
					break
				}
				prev = j
			}
		}

		if !curSuccess {
			prevSuccess = false
			continue
		}

		curEigenvalues, ok := extractEigenvalues(&squareMatrix)

		if !ok {
			continue
		}

		if prevSuccess && len(prevEigenvalues) == len(curEigenvalues) {
			success := true
			for i := 0; i < len(prevEigenvalues); i++ {
				if cmplx.Abs(prevEigenvalues[i]-curEigenvalues[i]) > STOP_THRESHOLD {
					//success = false
					break
				}
			}

			// additionally check elements under diagonal !
			if success {
				break
			}
		}

		prevEigenvalues = curEigenvalues
		prevSuccess = true
	}

	eigenvalues, _ := extractEigenvalues(&squareMatrix)
	return eigenvalues, extractEigenvectors(&squareMatrix, transform), iterations
}
