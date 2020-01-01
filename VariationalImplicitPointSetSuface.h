#ifndef VARIATIONALIMPLICITPOINTSETSURFACE_H
#define VARIATIONALIMPLICITPOINTSETSURFACE_H

#include <internal/STLContainerTraits.h>
#include <utils/NLoptimizer.h>
#include <Eigen/Dense>
#include <cmath>

template <typename Caller>
struct CubicRadialKernel
{
    typedef typename Caller::Scalar Scalar;
    CubicRadialKernel(const Caller* caller): mCaller(caller) {}
    Scalar operator()(const Eigen::Index& i, const Eigen::Index& j)
    {
        return std::pow((mCaller->mV.row(i) - mCaller->mV.row(j)).squaredNorm(), 1.5);
    }

    template <typename Derived>
    void D1(const Eigen::Index& i, const Eigen::Index& j, Eigen::MapBase<Derived>& derivate) const
    {
        auto p = (mCaller->mV.row(i) - mCaller->mV.row(j)).eval();
        derivate = p * (p.norm() * 3.0);
    }
    template <typename Derived>
    void D2(const Eigen::Index& i, const Eigen::Index& j, Eigen::MapBase<Derived>& hessian) const
    {
        auto p = (mCaller->mV.row(i) - mCaller->mV.row(j)).eval();
        auto len = p.norm();

        if(len < std::numeric_limits<typename Derived::Scalar>::epsilon())
        {
            hessian.setZero();
            return;
        }

        for(int l = 0; l < 3; l++)
        {
            for(int m = 0; m < 3; m++)
            {
                if(l == m)
                    hessian(l, m) = (p[l] * p[m] / len + len) * 3.0;
                else
                    hessian(l, m) = p[l] * p[m] / len * 3.0;
            }
        }

    }

    const Caller* mCaller;
};

template <typename Point3, template <typename> class KernelT>
class VariationalImplicitPointSetSurface
{
    typedef VariationalImplicitPointSetSurface<Point3, KernelT> Self;
    typedef KernelT<VariationalImplicitPointSetSurface<Point3, KernelT> > Kernel;
    friend Kernel;
public:
    typedef typename Point3::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    typedef Eigen::Matrix<Scalar, 1, 3> RowPoint3;

    template <typename Alloc>
    VariationalImplicitPointSetSurface(const std::vector<Point3, Alloc>& points, const Scalar& lambda);

    void Run();

public:
    template <typename DerivedM>
    void BuildM(Eigen::PlainObjectBase<DerivedM>& M);

    template <typename DerivedH, typename DerivedN>
    bool ComputeInitialNormals(const MatrixXd& J00, const MatrixXd& J10, const MatrixXd& J11, 
            const Scalar& lambda, Eigen::PlainObjectBase<DerivedH>& H, Eigen::PlainObjectBase<DerivedN>& N);

    template <typename DerivedN>
    bool ComputeNormals(Eigen::PlainObjectBase<DerivedN>& N, Scalar& energy);

    MatrixXd mV;
    MatrixXd mH;
    RowPoint3 mCenter;
    Scalar mScale;
    Scalar mLambda;
    Kernel mKernel;
};

template <typename Surface>
typename Surface::Scalar EnergyFuc( const std::vector<typename Surface::Scalar>& x, std::vector<typename Surface::Scalar>& grad, void* data)
{
    typedef typename Surface::Scalar Scalar;
    auto surface = reinterpret_cast<Surface*>(data);
    const Eigen::Index n = x.size() >> 1;

    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 2, Eigen::RowMajor>, 0, Eigen::OuterStride<2> > ex(&x[0], n, 2, Eigen::OuterStride<2>());    //eigen-style x
    Eigen::Matrix<Scalar, Eigen::Dynamic, 4> csMat(n, 4);     //cos and sin mat;
    csMat.col(0) = ex.col(0).array().cos().matrix();
    csMat.col(1) = ex.col(0).array().sin().matrix();
    csMat.col(2) = ex.col(1).array().cos().matrix();
    csMat.col(3) = ex.col(1).array().sin().matrix();

    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> g(n * 3);
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 3>, 0, Eigen::OuterStride<-1> > mg(g.data(), n, 3, Eigen::OuterStride<-1>(n));
    mg.col(0) = (csMat.col(0).array() * (csMat.col(2).array())).matrix();
    mg.col(1) = (csMat.col(0).array() * (csMat.col(3).array())).matrix();
    mg.col(2) = csMat.col(1);

    auto Hg = (surface->mH * g).eval();
    Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 3>, 0, Eigen::OuterStride<-1> > mHg(Hg.data(), n, 3, Eigen::OuterStride<-1>(n));

    if(!grad.empty())
    {
        grad.resize(n * 2);
        Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, 2, Eigen::RowMajor>, 0, Eigen::OuterStride<2> > eGrad(&grad[0], n, 2, Eigen::OuterStride<2>());
        eGrad.col(0) = 
           (-csMat.col(1).array() * csMat.col(2).array() * mHg.col(0).array()
            -csMat.col(1).array() * csMat.col(3).array() * mHg.col(1).array()
            +csMat.col(0).array() * mHg.col(2).array()).matrix() * 2.0;
        eGrad.col(1) = 
           (-csMat.col(0).array() * csMat.col(3).array() * mHg.col(0).array()
            +csMat.col(0).array() * csMat.col(2).array() * mHg.col(1).array()).matrix() * 2.0;
    }

    return g.transpose() * Hg;
}

template <typename Point3, template <typename> class Kernel>
template <typename Alloc>
VariationalImplicitPointSetSurface<Point3, Kernel>::VariationalImplicitPointSetSurface(const std::vector<Point3, Alloc>& points, const Scalar& lambda)
    :mLambda(lambda), mKernel(this)
{
    mV.resize(points.size(), 3);
    for(unsigned i = 0; i < points.size(); i++)
    {
        mV(i, 0) = points[i][0];
        mV(i, 1) = points[i][1];
        mV(i, 2) = points[i][2];
    }

    auto minCorner = mV.colwise().minCoeff();
    auto maxCorner = mV.colwise().maxCoeff();
    mCenter = (minCorner + maxCorner) * 0.5;
    Scalar maxLen = (maxCorner - minCorner).maxCoeff();
    mScale = 2.0 / maxLen;

    mV.rowwise() -= mCenter;
    mV *= mScale;
}

template <typename Point3, template <typename> class Kernel>
void VariationalImplicitPointSetSurface<Point3, Kernel>::Run()
{
    const auto dn = mV.rows() * 4;      //(dim + 1) x nPoints
    MatrixXd M(dn, dn);
    BuildM(M);

    MatrixXd N = MatrixXd::Zero(dn, 4);
    N.topLeftCorner(mV.rows(), 1).setOnes();
    N.topRightCorner(mV.rows(), 3) = mV;
    for(int i = 0; i < 3; i++)
        Eigen::Map<Eigen::Matrix<Scalar, -1, 1>>(N.data() + (dn + mV.rows()) * (i + 1), mV.rows()).setConstant(-1.0);

    MatrixXd A = MatrixXd::Zero(dn + 4, dn + 4);
    A.topLeftCorner(dn, dn) = M;
    A.topRightCorner(dn, 4) = N;
    A.bottomLeftCorner(4, dn) = N.transpose();
    auto invA = A.inverse();

    MatrixXd J00 = invA.topLeftCorner(mV.rows(), mV.rows());
    MatrixXd J10 = invA.block(mV.rows(), 0, mV.rows() * 3, mV.rows()); 
    MatrixXd J11 = invA.block(mV.rows(), mV.rows(), mV.rows() * 3, mV.rows() * 3);

    std::vector<Scalar> lambdas{0.0, 0.001, 0.01, 0.1, 1.0};
    typename Internal::STLVectorTraits<MatrixXd>::Vector normal_mats;
    std::vector<Scalar> energy(lambdas.size(), std::numeric_limits<Scalar>::max());

    for(unsigned i = 0; i < lambdas.size(); i++)
    {
        MatrixXd H, normal_mat;
        bool success = ComputeInitialNormals(J00, J10, J11, mLambda + lambdas[i], H, normal_mat);
        if(!i) mH = H;

        normal_mats.emplace_back(normal_mat);
        if(!success) continue;

        if(!ComputeNormals(normal_mat, energy[i])) continue;

        std::cout << "e: " << energy[i] << "\n\n";
    }
}

template <typename Point3, template <typename> class Kernel>
template <typename DerivedA>
void VariationalImplicitPointSetSurface<Point3, Kernel>::BuildM(Eigen::PlainObjectBase<DerivedA>& M)
{
    const auto nPoints = mV.rows();
    const auto matLength = nPoints * 4;

    //M00
    for(Eigen::Index i = 0; i < nPoints; i++)
    {
        for(Eigen::Index j = i; j < nPoints; j++)
        {
            M(i, j) = M(j, i) = mKernel(i, j);
        }
    }

    //M01 and M10
    auto colStride = Eigen::InnerStride<-1>(matLength * nPoints);       //across col
    auto rowStride = Eigen::InnerStride<-1>(nPoints);                   //across row
    for(Eigen::Index i = 0; i < nPoints; i++)
    {
        for(Eigen::Index j = 0; j < nPoints; j++)
        {
            Eigen::Map<RowPoint3, 0, Eigen::InnerStride<-1> > G(M.data() + ((j + nPoints) * matLength + i), colStride);
            mKernel.template D1<decltype(G)>(i, j, G);
            Eigen::Map<RowPoint3, 0, Eigen::InnerStride<-1> >(M.data() + (i * matLength + j + nPoints), rowStride) = G; 
        }
    }

    //M11
    auto stride = Eigen::Stride<-1, -1>(matLength * nPoints, nPoints);
    typedef Eigen::Matrix<Scalar, 3, 3> Matrix3d;
    for(Eigen::Index i = 0; i < nPoints; i++)
    {
        for(Eigen::Index j = i; j < nPoints; j++)
        {
            Eigen::Map<Matrix3d, 0, Eigen::Stride<-1, -1> > H(M.data() + ((nPoints + j) * matLength + nPoints + i), stride);
            mKernel.template D2<decltype(H)>(i, j, H);
            H *= -1.0;
            Eigen::Map<Matrix3d, 0, Eigen::Stride<-1, -1> >(M.data() + ((nPoints + i) * matLength + nPoints + j), stride) = H.transpose();
        }
    }
}

template <typename Point3, template <typename> class Kernel>
template <typename DerivedH, typename DerivedN>
inline bool VariationalImplicitPointSetSurface<Point3, Kernel>::ComputeInitialNormals(
        const MatrixXd& J00, const MatrixXd& J10,const MatrixXd& J11, const Scalar& lambda, Eigen::PlainObjectBase<DerivedH>& H, Eigen::PlainObjectBase<DerivedN>& N)
{
    if(lambda > std::numeric_limits<Scalar>::epsilon())
    {
        H = J11 - J10 * (DerivedN::Identity(mV.rows(), mV.rows()) + J00 * lambda).inverse() * J10.transpose() * lambda;
    }
    else
    {
        H = J11;
    }

    Eigen::SelfAdjointEigenSolver<DerivedN> solver(H);
    if(solver.info() != Eigen::Success) return false;

    auto eigen_vectors = solver.eigenvectors();
    N = Eigen::Map<DerivedN, 0, Eigen::Stride<0, 1> >(eigen_vectors.data(), mV.rows(), 3, Eigen::Stride<0, 1>());
    N.rowwise().normalize();
    return true;
}

template <typename Point3, template <typename> class Kernel>
template <typename DerivedN>
bool VariationalImplicitPointSetSurface<Point3, Kernel>::ComputeNormals(Eigen::PlainObjectBase<DerivedN>& N, Scalar& energy)
{
    std::vector<Scalar> x(mV.rows() * 2);
    for(Eigen::Index i = 0; i < mV.rows(); i++)
    {
        x[i*2]      = std::asin(N(i, 2));               //phi
        x[i*2 + 1]  = std::atan2(N(i, 1), N(i, 0));     //theta
    }

    std::vector<Scalar> lower(x.size(), -M_PI * 200.0);
    std::vector<Scalar> upper(x.size(), M_PI * 200.0);

    std::vector<Scalar> grad;
    grad.resize(x.size());

    std::cout << "enery:" << EnergyFuc<Self>(x, grad, this) << "\n";

    NLoptimizer optimizer;
    auto func = EnergyFuc<Self>;
    /*return optimizer(lower, upper, func, this, 1e-7, 100, x, energy) >= nlopt::SUCCESS;*/
    std::cout << "return: " << optimizer(lower, upper, func, this, 1e-7, 2000, x, energy) << std::endl;
    return true;
}
#endif
