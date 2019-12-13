package main

import (
	"cma-lab-go/matrix"
	"cma-lab-go/utils"
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"time"
)

func main() {
	fmt.Println("Cma lab on go! Here we go!")
	runtime.GOMAXPROCS(20) // number of CPUs [configure it later]
	rand.Seed(time.Now().UnixNano())

	// remember: const in go has another meaning - evaluated at compile time
	const SIZE = 3000 // untyped const

	//fmt.Printf("%.2f", utils.GenerateSquareMatrix(SIZE, -1, 1))

	lhs := utils.GenerateSquareMatrix(SIZE, -math.MaxInt64, math.MaxFloat64)
	rhs := utils.GenerateSquareMatrix(SIZE, -math.MaxInt64, math.MaxFloat64)

	start := time.Now()
	_, _ = matrix.MultiplyMatrices(lhs, rhs)
	fmt.Printf("SIZE = %v => time = %v milliseconds", SIZE,
		time.Since(start).Milliseconds())
}
