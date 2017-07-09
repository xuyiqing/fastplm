#include <vector>

#include "GPSolver.h"

typedef std::vector<arma::uword> OffsetT;

template <bool IsBalanced>
struct BalanceManager {
    GPSolver *ptr;
    
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
};

const bool kOffsetByTime  = true;
const bool kOffsetByIndiv = false;

template<> arma::vec
BalanceManager<false>::deflateVec(const arma::vec& vec, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::vec vec_(offset.size());
    for (int i = 0; i < offset.size(); i ++)
        vec_(i) = vec(offset[i]);
    return vec_;
}

template<> arma::vec
BalanceManager<true>::deflateVec(const arma::vec& vec, int id, bool byTime) {
    return vec;
}

template<> arma::mat
BalanceManager<false>::deflateMat(const arma::mat& mat, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::mat mat_(mat.n_rows, offset.size());
    for (int i = 0; i < mat.n_rows; i ++)
        mat_.row(i) = deflateVec(mat.row(i), id, byTime);
    return mat_;
}

template<> arma::mat
BalanceManager<true>::deflateMat(const arma::mat& mat, int id, bool byTime) {
    return mat;
}

template<> arma::vec
BalanceManager<false>::inflateVec(const arma::vec& vec, int id, bool byTime) {
    const auto& offset = byTime ? timeOffsets[id] : indivOffsets[id];
    arma::vec vec_ = arma::zeros(offset[offset.size() - 1] + 1);
    for (int i = 0; i < offset.size(); i ++)
        vec_(offset[i]) = vec(i);
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
BalanceManager<false>::BalanceManager(GPSolver *ptr_)
: ptr(ptr_), Ptoises({}), Piotses({}), timeOffsets({}), indivOffsets({})
{
    for (int i = 0; i < ptr->Y.n_rows; i ++) {
        OffsetT offset;
        for (int j = 0; j < ptr->Y.n_cols; j ++)
            if (!isnan(ptr->Y(i, j)))
                offset.push_back(j);
        
        timeOffsets.push_back(std::move(offset));
        Ptoises.push_back(makePtois(deflateMat(ptr->tois, i, true)));
    }
    
    for (int i = 0; i < ptr->Y.n_cols; i ++) {
        OffsetT offset;
        for (int j = 0; j < ptr->Y.n_rows; j ++)
            if (!isnan(ptr->Y(j, i)))
                offset.push_back(j);
        
        indivOffsets.push_back(std::move(offset));
        Piotses.push_back(makePtois(deflateMat(ptr->iots, i, false)));
    }
}

template<>
BalanceManager<true>::BalanceManager(GPSolver *ptr_): ptr(ptr_) {
    defaultPtois = makePtois(ptr->tois);
    defaultPiots = makePiots(ptr->iots);
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
    return Ptoises[time];
}

template<> const arma::mat&
BalanceManager<true>::getPiots(int _) {
    return defaultPiots;
}

template<bool IsBalanced>
void BalanceManager<IsBalanced>::MAP() {
    arma::cube data = arma::join_slices(ptr->Y, ptr->X);
    arma::cube copy;
    
    do {
        copy = data;
        
        for (int i = 0; i < data.n_slices; i ++) {
            arma::mat& panel = data.slice(i);
            
            for (int j = 0; j < ptr->indivCount; j ++) {
                auto col = deflateVec(panel.col(j), j, kOffsetByIndiv);
                panel.col(j) -= inflateVec(getPtois(j) * col, j, kOffsetByIndiv);
            }
            
            for (int j = 0; j < ptr->timeCount; j ++) {
                auto col = deflateVec(panel.row(j).t(), j, kOffsetByTime);
                panel.row(j) -= inflateVec(getPiots(j) * col, j, kOffsetByTime).t();
            }
        }
    } while (arma::accu(abs(data - copy)) > 1e-5);
    
    ptr->Y = data.slice(0);
    ptr->X = data.slices(1, data.n_slices - 1);
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
