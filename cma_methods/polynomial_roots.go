package cma_methods

import (
	"math"
)

const NEWTON_ACCURACY = 1e-12

func diff(polynomial []float64) []float64 {
	result := make([]float64, len(polynomial) - 1)
	for i := 1; i < len(polynomial); i++ {
		result[i-1] = polynomial[i] * float64(i)
	}
	return result
}

func value(polynomial []float64, point float64) float64 {
	var sum float64 = polynomial[0]
	power := point
	for i := 1; i < len(polynomial); i++ {
		sum += polynomial[i] * power
		power *= point
	}
	return sum
}

func checkFourier(polynomial []float64, point float64) bool {
	return value(diff(diff(polynomial)), point)*value(polynomial, point) > 0
}

func boundRoots(polynomial []float64) (float64, float64) {
	var max float64 = 0
	for i := 1; i < len(polynomial); i++ {
		max = math.Max(max, math.Abs(polynomial[i]))
	}
	high := math.Abs(polynomial[len(polynomial) - 1])
	return -1 - max / high, 1 + max / high
}

func bisection(polynomial []float64, left, right float64) (float64, float64) {
	rightSign := math.Signbit(value(polynomial, right))
	middle := (left + right) / 2

	if math.Signbit(value(polynomial, middle)) == rightSign {
		return left, middle
	}
	return middle, right
}

func newton(polynomial, differential []float64, point float64) float64 {
	return point - value(polynomial, point) / value(differential, point)
}

func appendRoot(roots, polynomial []float64, left, right float64) []float64 {
	if left == right {
		return roots
	}

	leftValue, rightValue := value(polynomial, left), value(polynomial, right)

	if leftValue == 0 {
		return append(roots, left)
	}

	if rightValue == 0 {
		return append(roots, right)
	}

	if math.Signbit(leftValue) == math.Signbit(rightValue) {
		return roots
	}

	middle := (left + right) / 2

	// do bisections until Fourier condition is true
	for ;!checkFourier(polynomial, middle); {
		left, right = bisection(polynomial, left, right)
		middle = (left + right) / 2
	}

	differential := diff(polynomial)
	prev := middle
	cur := newton(polynomial, differential, prev)
	for math.Abs(prev - cur) > NEWTON_ACCURACY {
		prev = cur
		cur = newton(polynomial, differential, prev)
	}

	return append(roots, cur)
}

func FindPolynomialRoots(polynomial []float64) []float64 {
	// linear polynomial
	deg := len(polynomial) - 1
	if deg == 1 {
		return []float64{-polynomial[0] / polynomial[1]}
	}

	var roots []float64
	leftBound, rightBound := boundRoots(polynomial)
	extremes := FindPolynomialRoots(diff(polynomial))

	if len(extremes) == 0 {
		return appendRoot(roots, polynomial, leftBound, rightBound)
	}

	roots = appendRoot(roots, polynomial, leftBound, extremes[0])
	for i := 1; i < len(extremes); i++ {
		roots = appendRoot(roots, polynomial, extremes[i - 1], extremes[i])
	}
	roots = appendRoot(roots, polynomial, extremes[len(extremes) - 1], rightBound)

	return roots
}
