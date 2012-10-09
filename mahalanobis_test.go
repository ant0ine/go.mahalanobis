package mahalanobis

import (
	//    "fmt"
	"github.com/skelterjohn/go.matrix"
	"math"
	"testing"
)

func TestMeanVector(t *testing.T) {

	points, err := matrix.ParseMatlab("[1 1 1;1 1 1]")
	if err != nil {
		t.Fatal(err)
	}

	expected, err := matrix.ParseMatlab("[1;1]")
	if err != nil {
		t.Fatal(err)
	}

	result := MeanVector(points)
	if !matrix.Equals(result, expected) {
		t.Error()
	}

	points, err = matrix.ParseMatlab("[0 1 2;0 2 4]")
	if err != nil {
		t.Fatal(err)
	}

	expected, err = matrix.ParseMatlab("[1;2]")
	if err != nil {
		t.Fatal(err)
	}

	result = MeanVector(points)
	if !matrix.Equals(result, expected) {
		t.Error()
	}
}

func TestCovarianceMatrix(t *testing.T) {

	// no (co)variance
	// R: var(cbind(c(1, 1), c(1, 1)))
	points, err := matrix.ParseMatlab("[1 1;1 1]")
	if err != nil {
		t.Fatal(err)
	}

	expected, err := matrix.ParseMatlab("[0 0;0 0]")
	if err != nil {
		t.Fatal(err)
	}

	result := CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.Equals(result, expected) {
		t.Error()
	}

	// diagonale case
	// R: var(cbind(c(0, 4, 2, 2), c(2, 2, 0, 4)))
	points, err = matrix.ParseMatlab("[0 4 2 2;2 2 0 4]")
	if err != nil {
		t.Fatal(err)
	}

	// R: var(cbind(c(0, 4, 2, 2), c(2, 2, 0, 4)))
	expected = matrix.MakeDenseMatrix([]float64{2.66, 0, 0, 2.66}, 2, 2)

	result = CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.ApproxEquals(result, expected, 0.01) {
		t.Error()
	}

	// another case
	// R: var(cbind(c(9, 3, 5), c(3, 4, 1)))
	points, err = matrix.ParseMatlab("[9 3 5;3 4 1]")
	if err != nil {
		t.Fatal(err)
	}

	expected = matrix.MakeDenseMatrix([]float64{9.33, -0.66, -0.66, 2.33}, 2, 2)

	result = CovarianceMatrix(points)
	//fmt.Println("covariance:\n", result)
	if !matrix.ApproxEquals(result, expected, 0.01) {
		t.Error()
	}
}

func TestDistance(t *testing.T) {

	// R:
	// x = cbind(c(9, 3, 5), c(3, 4, 1))
	// mahalanobis(c(1,1), colMeans(x), var(x))
	points, err := matrix.ParseMatlab("[9 3 5;3 4 1]")
	if err != nil {
		t.Fatal(err)
	}

	target, err := matrix.ParseMatlab("[1;1]")
	if err != nil {
		t.Fatal(err)
	}

	square, err := DistanceSquare(points, target)
	if err != nil {
		t.Fatal(err)
	}

	if math.Abs(square-4.08) > 0.01 {
		t.Error()
	}
}
