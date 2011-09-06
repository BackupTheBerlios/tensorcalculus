/*
 * Copyright (C) 2011 Stefan Handschuh
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "TensorTrainRepresentation.hpp" // includes TensorRepresentation.hpp, TensorChainRepresentation.hpp, LapackInterface2.hpp and BlasInterface.hpp




namespace TensorCalculus {

/**
 * This class handles the computation of the gradient of an
 * arbitrary tensor representation where the objective function is
 *     \f$f_a(x) = \|a-x\|\f$
 * This results in a slow implementation in most cases
 * (except CP-tensor for instance) such that a separate
 * implementation for tensor networks is highly recommended.
 * When doing this, there is no need to derive from this class
 * although it is fine to do this.
 *
 * An other way of adding optimized implementation is to add
 * implemented templates.
 *
 */

template<typename T>
class DistanceFunctionGradient;

template<typename T>
class DistanceFunctionGradient< TensorRepresentation<T> > {

protected:
    TensorRepresentation<T> a; // fixed representation a

    TensorRepresentation<T> x; // approximated representation

    int d;

    T normA;

    void init(const TensorRepresentation<T> &a, const TensorRepresentation<T> &x) {
    	this->a = a;
    	this->x = x;
    	this->d = a.getD();
    	normA = l2norm(a);
    }

public:
    DistanceFunctionGradient() {
    	this->d = 0;
    	this->normA = 0;
    }

    DistanceFunctionGradient(const TensorRepresentation<T> &a, const TensorRepresentation<T> &x) {
    	init(a, x);
    }

    /**
     * Computes the gradient at $x and creates the result which is returned.
     */
    std::vector< std::vector<T> > operator () (const std::vector< std::vector<T> >& x) {
        std::vector< std::vector<T> > result(d);

        for (int n = 0; n < d; n++) {
           result[n] = std::vector<T>(x[n].size());
        }
        operator () (x, result);
        return result;
    }

    /**
     * Evaluates the gradient at $x and stores the result in $result.
     */
    void operator () (const std::vector< std::vector<T> >& x, std::vector< std::vector<T> >& result) {
        this->x.setV(x); // since we are using the methods of this->x

        std::vector< std::vector<T> > scalarComponentsXX; // components <x_i(j_1), x_i(j_2)>
        std::vector< std::vector<T> > scalarComponentsAX; // components <a_i(j_1), x_i(j_2)>

        this->x.tensorScalarProductComponents(this->x, scalarComponentsXX); // compute components
        a.tensorScalarProductComponents(this->x, scalarComponentsAX);       // compute components

        for (int n = 0; n < d; n++) {
        	using namespace VectorOperators;
        	result[n] *= 0.0; // reset the result

            int dimension = this->x.getComponentDimension(n); // current component dimension

            std::vector<int> partialSummation = this->x.getPartialSummations(n);

            Index index2(partialSummation);

            std::vector<int> partialSummationB = a.getPartialSummations(n);

            Index index3(partialSummationB);

            for (Index index1(partialSummation); !index1.end(); ++index1) {
				for (index2.begin(); !index2.end(); ++index2) {
					const T factor = this->x.tensorScalarProduct(scalarComponentsXX, n, index2, index1, this->x);

					Blas<T>::axpy(dimension, factor, &(this->x.getV()[n][index2*dimension]), 1, &result[n][index1*dimension], 1);
				}
				for (index3.begin(); !index3.end(); ++index3) {
				    const T factor = -a.tensorScalarProduct(scalarComponentsAX, n, index3, index1, this->x);

				    Blas<T>::axpy(dimension, factor, &a.getV()[n][index3*dimension], 1, &result[n][index1*dimension], 1);
				}
            }
        }
    }

}; // end of class


template<typename T>
class DistanceFunctionGradient< TensorChainRepresentation<T> > {
private:
	TensorChainRepresentation<T> a;

	TensorChainRepresentation<T> x;

	int d;

	T normA;

	std::vector< std::vector<T> > ascA;

	std::vector< std::vector<T> > descA;

	std::vector< std::vector<T> > ascB;

	std::vector< std::vector<T> > descB;

	std::vector<T> dummyA;

	std::vector<T> dummyB;

	std::vector<T> dummy;

protected:
	void init(const TensorChainRepresentation<T> &a, const TensorChainRepresentation<T> &x) {
		this->a = a;
		this->x = x;
		this->d = a.getD();
		this->normA = l2norm(a);
		ascA.resize(this->d-1);
		descA.resize(this->d-1);
		ascB.resize(this->d-1);
		descB.resize(this->d-1);
	}

public:
    DistanceFunctionGradient() {
    	this->d = 0;
    	this->normA = 0;
    }

    DistanceFunctionGradient(const TensorChainRepresentation<T> &a, const TensorChainRepresentation<T> &x) {
        this->a = a;
        this->x = x;
        this->d = a.getD();
        this->normA = l2norm(a);
    }

    /**
     * Computes the gradient at $x and creates the result which is returned.
     *
     * We only need to override this due to cpp compiler complaints.
     */
    std::vector< std::vector<T> > operator () (const std::vector< std::vector<T> >& x) {
        std::vector< std::vector<T> > result(this->d);

        for (int n = 0; n < this->d; n++) {
            result[n] = std::vector<T>(x[n].size());
        }
        operator () (x, result);
        return result;
    }

    /**
     * Evaluates the gradient at $x and stores the result in $result.
     */
    void operator () (const std::vector< std::vector<T> >& x, std::vector< std::vector<T> >& result) {
    	this->x.setV(x); // since we are using the methods of this->x

    	int m; // rows

    	int n; // cols

    	int k; // rows

    	dummy.resize(std::max(this->x.getSubsequentMaxRankProduct(this->x), this->x.getSubsequentMaxRankProduct(this->a)));

    	/* Compute the desc structure */
    	this->x.fill_A_mu(this->a, 0, descA[0]);
    	this->x.fill_B_mu(0, descB[0]);
    	for (int mu = 1; mu < this->d-1; mu++) {
    		const int mu_plus_1 = (mu + 1) % this->d;

    		this->x.fill_A_mu(this->a, mu, dummyA);
    		this->x.fill_B_mu(mu, dummyB);

    		m = this->a.getSummation(0)*this->x.getSummation(0);
    		n = this->a.getSummation(mu_plus_1)*this->x.getSummation(mu_plus_1);
    		k = this->a.getSummation(mu)*this->x.getSummation(mu);
    		descA[mu].resize(m*n);

    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &descA[mu-1][0], &dummyA[0], 0.0, &descA[mu][0]); // descA[mu] = descA[mu-1]*dummyA

    		m = sqr(this->x.getSummation(0));
    		n = sqr(this->x.getSummation(mu_plus_1));
    		k = sqr(this->x.getSummation(mu));

    		descB[mu].resize(m*n);
    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &descB[mu-1][0], &dummyB[0], 0.0, &descB[mu][0]); // descB[mu] = descB[mu-1]*dummyB;
    	}

    	/* Compute the asc structure */
    	this->x.fill_A_mu(this->a, this->d-1, ascA[this->d-2]);
    	this->x.fill_B_mu(this->d-1, ascB[this->d-2]);
    	for (int mu = this->d-2; mu > 0; mu--) {
    		this->x.fill_A_mu(this->a, mu, dummyA);
    		this->x.fill_B_mu(mu, dummyB);

    		m = this->a.getSummation(mu)*this->x.getSummation(mu);
    		n = this->a.getSummation(0)*this->x.getSummation(0);
    		k = this->a.getSummation(mu+1)*this->x.getSummation(mu+1);

    		ascA[mu-1].resize(m*n);
    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &dummyA[0], &ascA[mu][0], 0.0, &ascA[mu-1][0]); // ascA[mu-1] = dummyA*ascA[mu];

    		m = sqr(this->x.getSummation(mu));
    		n = sqr(this->x.getSummation(0));
    		k = sqr(this->x.getSummation(mu+1));
    		ascB[mu-1].resize(m*n);
    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &dummyB[0], &ascB[mu][0], 0.0, &ascB[mu-1][0]); // descB[mu-1] = dummyB*ascB[mu];
    	}

    	/* Compute the gradient */
    	int dimension;

    	/* Special case for the first */
    	dimension = this->x.getComponentDimension(0);
    	this->x.reshape_into(this->a, 0, 1, ascA[0], dummyA); //dummyA = reshape(ascA[0]));
    	this->x.reshape_into(0, 1, ascB[0], dummyB); //dummyB = reshape(ascB[0]));

    	n = this->a.getSummation(0)*this->a.getSummation(1);
    	k = this->x.getSummation(0)*this->x.getSummation(1);

    	Blas<T>::gemm('N', 'T', dimension, n, k, -1.0, &(this->a.getV()[0][0]), &dummyA[0], 0.0, &result[0][0]); // result[0] = - a[0]*ascA[0]^T
    	Blas<T>::gemm('N', 'T', dimension, k, k, 1.0, &x[0][0], &dummyB[0], 1.0, &result[0][0]); // result[0] += x[0]*ascB[0]^T;

    	/* Standard cases */
    	for (int mu = 1; mu < this->d-1; mu++) {
    		dimension = this->x.getComponentDimension(mu);

    		m = this->a.getSummation(mu)*this->x.getSummation(mu);
    		n = this->a.getSummation(mu+1)*this->x.getSummation(mu+1);
    		k = this->a.getSummation(0)*this->x.getSummation(0);

    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &ascA[mu][0], &descA[mu-1][0], 0.0, &dummy[0]); // dummy = ascA[mu-1]*descA[mu]
    		this->x.reshape_into(this->a, mu, mu+1, dummy, dummyA);

    		m = sqr(this->x.getSummation(mu));
    		n = sqr(this->x.getSummation(mu+1));
    		k = sqr(this->x.getSummation(0));

    		Blas<T>::gemm('N', 'N', m, n, k, 1.0, &ascB[mu][0], &descB[mu-1][0], 0.0, &dummy[0]); // dummy = ascB[mu-1]*descB[mu]
    		this->x.reshape_into(mu, mu+1, dummy, dummyB);

    		n = this->a.getSummation(mu)*this->a.getSummation(mu+1);
    		k = this->x.getSummation(mu)*this->x.getSummation(mu+1);
    		Blas<T>::gemm('N', 'T', dimension, n, k, -1.0, &(this->a.getV()[mu][0]), &dummyA[0], 0.0, &result[mu][0]); // result[mu] = -a[mu]*dummyA^T
    		Blas<T>::gemm('N', 'T', dimension, k, k, 1.0, &x[mu][0], &dummyB[0], 1.0, &result[mu][0]); // result[mu] += x[mu]*dummyB^T;
    	}

    	/* Special case for the last */
    	dimension = this->x.getComponentDimension(this->d-1);
    	this->x.reshape_into(this->a, this->d-1, 0, descA[this->d-2], dummyA); //dummyA = reshape(ascA[0]));
    	this->x.reshape_into(this->d-1, 0, descB[this->d-2], dummyB); //dummyB = reshape(ascB[0]));

    	n = this->a.getSummation(this->d-1)*this->a.getSummation(0);
    	k = this->x.getSummation(this->d-1)*this->x.getSummation(0);

    	Blas<T>::gemm('N', 'T', dimension, n, k, -1.0, &(this->a.getV()[this->d-1][0]), &dummyA[0], 0.0, &result[this->d-1][0]); // result[this->d-1] = -a[this->d-1]*descA[this->d-2]^T
    	Blas<T>::gemm('N', 'T', dimension, k, k, 1.0, &x[this->d-1][0], &dummyB[0], 1.0, &result[this->d-1][0]); // result[this->d-1] += x[this->d-1]*descB[this->d-2]^T;
    }
};


/**
 * We do not need to implement this class since is is already implemented with the tensor chain
 * format where we have kept in mind that the last edge is omitted in tensor trains.
 *
 * This class is only there for cosmetic reasons such that it can be directly used instead of
 * the TensorChain.
 */
template<typename T>
class DistanceFunctionGradient< TensorTrainRepresentation<T> > :
public DistanceFunctionGradient< TensorChainRepresentation<T> > {

public:
	DistanceFunctionGradient() {
		this->d = 0;
		this->normA = 0;
	}

	DistanceFunctionGradient(const TensorTrainRepresentation<T> &a, const TensorTrainRepresentation<T> &x) {
	    init(a, x);
	}

	using DistanceFunctionGradient< TensorChainRepresentation<T> >::operator();
};

} // end of namespace
