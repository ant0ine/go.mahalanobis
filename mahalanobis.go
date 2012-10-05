// Naive implementation of the Mahalanobis distance using go.matrix
// (https://en.wikipedia.org/wiki/Mahalanobis_distance)
//
// This is me learning Go, it's probably broken, don't use it.
//
// Example:
//
//      package main
//
//      import (
//          "fmt"
//          "github.com/skelterjohn/go.matrix"
//          "github.com/ant0ine/go.mahalanobis"
//      )
//
//      func main() {
//
//          points, err := matrix.ParseMatlab("[1 4 3 4;4 2 3 4]")
//          if err != nil {
//          panic(err)
//          }
//          fmt.Println("4 points:\n", points)
//
//          target, err := matrix.ParseMatlab("[3;4]")
//          if err != nil {
//              panic(err)
//          }
//          fmt.Println("the target point:\n", target)
//
//          distance, err := mahalanobis.Distance(points, target)
//          if err != nil {
//              panic(err)
//          }
//          fmt.Println("Mahalanobis distance=", distance)
//      }
package mahalanobis

import (
//    "fmt"
    "math"
    "github.com/skelterjohn/go.matrix"
)

// Given a set a points, return the mean vector.
func MeanVector(points *matrix.DenseMatrix) *matrix.DenseMatrix {
    mean := matrix.Zeros(points.Rows(), 1)
    for i := 0; i < points.Rows(); i++ {
        sum := 0.0
        for j := 0; j < points.Cols(); j++ {
            sum += points.Get(i, j)
        }
        mean.Set(i, 0, sum / float64(points.Cols()))
    }
    return mean
}

func sample_covariance_matrix(points, mean *matrix.DenseMatrix) *matrix.DenseMatrix {
    dim := points.Rows()
    cov := matrix.Zeros(dim, dim)
    for i := 0; i < dim; i++ {
        for j := 0; j < dim; j++ {
            if i > j {
                // symetric matrix
                continue
            }
            // TODO in go routines ?
            sum := 0.0
            for k := 0; k < points.Cols(); k++ {
                sum += (points.Get(i, k) - mean.Get(i, 0)) * (points.Get(j, k) - mean.Get(j, 0))
            }

            // this is the sample covariance, divide by (N - 1)
            covariance := sum / ( float64(points.Cols() - 1))

            cov.Set(i, j, covariance)
            // symetric matrix
            cov.Set(j, i, covariance)

        }
    }
    return cov
}

// Return the covariance matrix for this set of points (sample covariance is used)
func CovarianceMatrix(points *matrix.DenseMatrix) *matrix.DenseMatrix {
    mean := MeanVector(points)
    return sample_covariance_matrix(points, mean)
}

// Return the square of the Mahalanobis distance
func DistanceSquare(points, target *matrix.DenseMatrix) (float64, error) {

    // TODO check the dimensions

    mean := MeanVector(points)
    //fmt.Println("mean:\n", mean)

    delta := target.Copy()
    delta.SubtractDense(mean)
    //fmt.Println("delta:\n", delta)

    cov := sample_covariance_matrix(points, mean)
    //fmt.Println("covariance:\n", cov)

    inv, err := cov.Inverse()
    if err != nil {
        return 0, err // XXX
    }
    //fmt.Println("inverse covariance:\n", inv)

    product1, err := inv.TimesDense(delta)
    if err != nil {
        return 0, err // XXX
    }
    delta_t := delta.Transpose()
    product2, err := delta_t.TimesDense(product1)
    if err != nil {
        return 0, err // XXX
    }

    return product2.Get(0,0),  nil
}

// Return the Mahalanobis distance
func Distance(points, target *matrix.DenseMatrix) (float64, error) {
    square, err := DistanceSquare(points, target)
    if err != nil {
        return 0, err // XXX
    }
    return math.Sqrt(square), nil

}
