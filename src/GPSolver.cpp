#include <vector>

#include "GPSolver.h"

typedef std::vector<arma::uword> OffsetT;

template <typename T>
void deflateMem(const T *src, T *dest, const OffsetT &offset) {
    for (int i = 0; i < offset.size(); i ++)
        dest[i] = src[offset[i]];
}

template <typename T>
void inflateMem(const T *src, T *dest, const OffsetT &offset) {
    for (int i = 0; i < offset.size(); i ++)
        dest[offset[i]] = src[i];
}

template <bool IsBalanced>
struct BalanceManager {
    GPSolver *solver;
    
    arma::mat defaultPtois;
    arma::mat defaultPiots;
    
    std::vector<arma::mat> Ptoises;
    std::vector<arma::mat> Piotses;
    std::vector<OffsetT> timeOffsets;
    std::vector<OffsetT> indivOffsets;
    
public:
    arma::vec deflateVec(const arma::vec& vec, int id, bool byTime);
    arma::mat deflateMat(const arma::mat& mat, int id, bool byTime);
    arma::vec inflateVec(const arma::vec& vec, int id, bool byTime);
    
    BalanceManager(GPSolver *ptr);
    
    const arma::mat& getPtois(int indiv);
    const arma::mat& getPiots(int time);
    
    void MAP();
    
    void flatenUnbalancedMatrix(const double *src, double *dest);
    arma::vec flatenAndSolve();
};

const bool kOffsetByTime  = true;
const bool kOffsetByIndiv = false;

template<> arma::vec
BalanceManager<false>::deflateVec(const arma::vec& vec, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::vec vec_(offset.size());
    deflateMem(vec.memptr(), vec_.memptr(), offset);
    return vec_;
}

template<> arma::vec
BalanceManager<true>::deflateVec(const arma::vec& vec, int id, bool byTime) {
    return vec;
}

template<> arma::mat
BalanceManager<false>::deflateMat(const arma::mat& mat, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::mat mat_(offset.size(), mat.n_cols);
    for (int i = 0; i < mat.n_cols; i ++)
        mat_.col(i) = deflateVec(mat.col(i), id, byTime);
    
    return mat_;
}

template<> arma::mat
BalanceManager<true>::deflateMat(const arma::mat& mat, int id, bool byTime) {
    return mat;
}

template<> arma::vec
BalanceManager<false>::inflateVec(const arma::vec& vec, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::vec vec_(offset[offset.size() - 1] + 1);
    inflateMem(vec.memptr(), vec_.memptr(), offset);
    return vec_;
}

template<> arma::vec
BalanceManager<true>::inflateVec(const arma::vec& vec, int id, bool byTime) {
    return vec;
}

inline arma::mat makePtois(const arma::mat& tois) {
    return tois * arma::inv(tois.t() * tois) * tois.t();
}

inline arma::mat makePiots(const arma::mat& iots) {
    return makePtois(iots);
}

template<>
BalanceManager<false>::BalanceManager(GPSolver *solver_)
: solver(solver_), Ptoises({}), Piotses({}), timeOffsets({}), indivOffsets({})
{
    // time-indexed offsets for IOTS data
    for (int i = 0; i < solver->Y.n_rows; i ++) {
        OffsetT offset;
        for (int j = 0; j < solver->Y.n_cols; j ++)
            if (!isnan(solver->Y(i, j)))
                offset.push_back(j);
        
        timeOffsets.push_back(std::move(offset));
        
        auto mat_ = deflateMat(solver->iots, i, kOffsetByTime);
        Piotses.push_back(makePiots(mat_));
    }
    
    // indiv-indexed offsets for TOIS data
    for (int i = 0; i < solver->Y.n_cols; i ++) {
        OffsetT offset;
        for (int j = 0; j < solver->Y.n_rows; j ++)
            if (!isnan(solver->Y(j, i)))
                offset.push_back(j);
        
        indivOffsets.push_back(std::move(offset));
        
        auto mat_ = deflateMat(solver->tois, i, kOffsetByIndiv);
        Ptoises.push_back(makePtois(mat_));
    }
}

template<>
BalanceManager<true>::BalanceManager(GPSolver *ptr_): solver(ptr_) {
    defaultPtois = makePtois(solver->tois);
    defaultPiots = makePiots(solver->iots);
}

template<> const arma::mat&
BalanceManager<false>::getPtois(int indiv) {
    return Ptoises[indiv];
}
template<> const arma::mat&
BalanceManager<true>::getPtois(int _) {
    return defaultPtois;
}

template<> const arma::mat&
BalanceManager<false>::getPiots(int time) {
    return Piotses[time];
}

template<> const arma::mat&
BalanceManager<true>::getPiots(int _) {
    return defaultPiots;
}

template<bool IsBalanced>
void BalanceManager<IsBalanced>::MAP() {
    arma::cube data = arma::join_slices(solver->Y, solver->X);
    arma::cube copy;
    
    do {
        copy = data;
        
        for (int i = 0; i < data.n_slices; i ++) {
            arma::mat& panel = data.slice(i);
            
            for (int j = 0; j < solver->indivCount; j ++) {
                auto col1 = deflateVec(panel.col(j), j, kOffsetByIndiv);
                auto col2 = inflateVec(getPtois(j) * col1, j, kOffsetByIndiv);
                panel.col(j) -= col2;
            }
            
            for (int j = 0; j < solver->timeCount; j ++) {
                auto col1 = deflateVec(panel.row(j).t(), j, kOffsetByTime);
                auto col2 = inflateVec(getPiots(j) * col1, j, kOffsetByTime);
                panel.row(j) -= col2.t();
            }
        }
    } while (arma::accu(abs(data - copy)) > 1e-5);
    
    solver->Y = data.slice(0);
    solver->X = data.slices(1, data.n_slices - 1);
}

template<bool IsBalanced>
void BalanceManager<IsBalanced>::flatenUnbalancedMatrix(const double *src, double *dest) {
    for (const auto &indivOffset : indivOffsets) {
        deflateMem(src, dest, indivOffset);
        src += indivOffset.size();
        dest += solver->timeCount;
    }
}

template <> arma::vec
BalanceManager<false>::flatenAndSolve() {
    size_t obsCount = 0;
    for (const auto &xs : timeOffsets)
        obsCount += xs.size();
    
    std::shared_ptr<double> memory(new double[obsCount * (solver->paramCount + 1)],
                                   std::default_delete<double[]>());
    
    double *ptr = memory.get();
    flatenUnbalancedMatrix(solver->Y.memptr(), ptr);
    ptr += obsCount;
    
    for (int i = 0; i < solver->paramCount; i ++) {
        flatenUnbalancedMatrix(solver->X.slice(i).memptr(), ptr);
        ptr += obsCount;
    }
    
    arma::vec vecY(memory.get(),obsCount, 1, false);
    arma::mat matX(memory.get() + obsCount, obsCount, solver->paramCount, false);
    
    return solve(matX, vecY);
}

template <> arma::vec
BalanceManager<true>::flatenAndSolve() {
    auto obsCount = solver->timeCount * solver->indivCount;
    
    const arma::mat matX(solver->X.memptr(), obsCount, solver->paramCount, false);
    const arma::vec vecY(solver->Y.memptr(), obsCount, 1, false);
    
    return solve(matX, vecY);
}

GPSolver::GPSolver(arma::cube X, arma::mat Y,
                     arma::mat tois_, arma::mat iots_,
                     bool isBalanced, bool withFixedEffects):
    paramCount(X.n_slices), timeCount(X.n_rows), indivCount(X.n_cols),
    X(std::move(X)), Y(std::move(Y)),
    tois(std::move(tois_)), iots(std::move(iots_)),
    isBalanced(isBalanced)
{
    if (withFixedEffects) {
        tois.insert_cols(0, 1);
        tois.col(0) = arma::ones(timeCount);
        
        iots.insert_cols(0, 1);
        iots.col(0) = arma::ones(indivCount);
    }
}

void GPSolver::MAP() {
    if (isBalanced) {
        BalanceManager<true> manager(this);
        manager.MAP();
    }
    else {
        BalanceManager<false> manager(this);
        manager.MAP();
    }
}

arma::colvec GPSolver::compute() {
    MAP();
    
    const arma::mat matX(X.memptr(), timeCount * indivCount, paramCount, false);
    const arma::colvec vecY(Y.memptr(), timeCount * indivCount, 1, false);
    
    return solve(matX, vecY);
}
