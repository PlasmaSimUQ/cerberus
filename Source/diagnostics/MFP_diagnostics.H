#ifdef PYTHON
#ifndef MFP_DIAGNOSTICS_H
#define MFP_DIAGNOSTICS_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include "MFP_polyspline.H"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace amrex;

void calc_gradient(const Box& box,
                   FArrayBox& data,
                   FArrayBox& slope,
                   const int d);

std::string print_FAB(const FArrayBox& src, int n);

void plot_FAB_1d(const Box& box,
                 const FArrayBox& src,
                 const Vector<int> N,
                 std::string title,
                 bool block=false);

void plot_FAB_1d(const FArrayBox& src,
                 const Vector<int> N,
                 bool block,
                 std::string title);

void plot_FAB_1d(const FArrayBox& src,
                 bool block,
                 std::string title);

struct PlotData2D {
        Vector<float> data;
        int nr, nc;
        Vector<float> extents;
};

template<typename T>
PlotData2D get_2d_data(const Box& box,
                 const T& src,
                 const int n,
                 bool log=false)
{
    PlotData2D pd;

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real const> const& src4 = src.array();

    Vector<float> &data = pd.data;

    for   (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            float d = (float)src4(i,j,0,n);
            if (log) {
                if (std::abs(d) <= 0.0) {
                    d = std::nanf("1");
                }
                d = std::log10(std::abs(d));
            }
            data.push_back(d);
        }
    }

    pd.nr = hi.y-lo.y + 1;
    pd.nc = hi.x-lo.x + 1;

    pd.extents = {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f};

    return pd;
}

template<typename T>
void plot_FAB_2d(const Box &box,
                 const T& src,
                 const Vector<int> N,
                 std::string title="",
                 bool log=false,
                 bool block=false)
{
    plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    npr = (int)std::floor(std::sqrt(N.size()));
    npc = npr + N.size() - npr*npr;

    //        npc = std::floor(std::sqrt(N.size()));
    //        npr = npc + N.size() - npc*npc;

    for (int n=0; n<N.size(); ++n) {

        PlotData2D pd = get_2d_data(box, src, n, log);


        plt::subplot(npr,npc,n+1);
        plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);
        plt::colorbar();
        plt::title(std::to_string(N[n]));

    }

    plt::suptitle(title);


    if (block)
        plt::show(block);

    return;
}

template<typename T>
void plot_FAB_2d(const Box &box,
                 const T& src,
                 const int n,
                 std::string title="",
                 bool log=false,
                 bool block=false)
{
    plt::figure_size(1200, 780);

    PlotData2D pd = get_2d_data(box, src, n, log);

    plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);
    plt::colorbar();
    plt::title(title);

    if (block)
        plt::show(block);
    return;
}

template<typename T>
void plot_FAB_2d(const T& src,
                 const int n,
                 std::string title="",
                 bool log=false,
                 bool block=false)
{
    plt::figure_size(1200, 780);

    PlotData2D pd = get_2d_data(src.box(), src, n, log);

    plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);
    plt::colorbar();
    plt::title(title);

    if (block)
        plt::show(block);
    return;
}


template<typename T>
void plot_FAB_2d(const T& src,
                 const int n,
                 const Box &plot_box,
                 std::string title="",
                 bool log=false,
                 bool block=false)
{
    plt::figure_size(1200, 780);

    PlotData2D pd = get_2d_data(src.box(), src, n, log);

    plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);

    Vector<float> x(5), y(5);

    x[0] = plot_box.smallEnd()[0] - 0.5f;
    y[0] = plot_box.smallEnd()[1] - 0.5f;

    x[1] = plot_box.bigEnd()[0] + 0.5f;
    y[1] = plot_box.smallEnd()[1] - 0.5f;

    x[2] = plot_box.bigEnd()[0] + 0.5f;
    y[2] = plot_box.bigEnd()[1] + 0.5f;

    x[3] = plot_box.smallEnd()[0] - 0.5f;
    y[3] = plot_box.bigEnd()[1] + 0.5f;

    x[4] = plot_box.smallEnd()[0] - 0.5f;
    y[4] = plot_box.smallEnd()[1] - 0.5f;

    plt::plot(x, y,"k-");

    plt::colorbar();
    plt::title(title);

    if (block)
        plt::show(block);
    return;
}

template<typename T>
void plot_FAB_2d_vec(const T& src,
                 const Array<int,2> n,
                 const Box &plot_box,
                 std::string title="",
                 bool log=false,
                 bool block=false)
{
    plt::figure_size(1200, 780);

    PlotData2D pd1 = get_2d_data(src.box(), src, n[0], log);
    PlotData2D pd2 = get_2d_data(src.box(), src, n[1], log);

    std::vector<float> x(pd1.nr*pd1.nc), y(pd1.nr*pd1.nc);

    float start_x = pd1.extents[0] + 0.5f;
    float start_y = pd1.extents[2] + 0.5f;
    for (int j=0; j<pd1.nr; ++j) {
        for (int i=0; i<pd1.nc; ++i) {
            x[i+j*pd1.nc] = start_x + i;
            y[i+j*pd1.nc] = start_y + j;
//            x[j+i*pd1.nr] = start_x + i;
//            y[j+i*pd1.nr] = start_y + j;
        }
    }

    plt::imshow(pd1.data.data(), pd1.nr, pd1.nc, 1, {{"origin","lower"}}, pd1.extents);

    plt::quiver(x, y, pd1.data, pd2.data);

    x.resize(5);
    y.resize(5);

    x[0] = plot_box.smallEnd()[0] - 0.5f;
    y[0] = plot_box.smallEnd()[1] - 0.5f;

    x[1] = plot_box.bigEnd()[0] + 0.5f;
    y[1] = plot_box.smallEnd()[1] - 0.5f;

    x[2] = plot_box.bigEnd()[0] + 0.5f;
    y[2] = plot_box.bigEnd()[1] + 0.5f;

    x[3] = plot_box.smallEnd()[0] - 0.5f;
    y[3] = plot_box.bigEnd()[1] + 0.5f;

    x[4] = plot_box.smallEnd()[0] - 0.5f;
    y[4] = plot_box.smallEnd()[1] - 0.5f;

    plt::plot(x, y,"k-");

    plt::colorbar();
    plt::title(title);

    if (block)
        plt::show(block);
    return;
}

#ifdef AMREX_USE_EB
PlotData2D get_2d_data(const Box& box,
                 const EBCellFlagFab& src,
                 const int n,
                 bool log=false);

void plot_FAB_2d(const Box &box,
                 const EBCellFlagFab& src,
                 std::string title="",
                 bool block=false);

void plot_FAB_2d(const EBCellFlagFab& src,
                 std::string title="",
                 bool block=false);

void plot_FAB_eb(const FArrayBox& src,
                 const FArrayBox &pts,
                 const int n,
                 std::string title="",
                 bool block=false);

void plot_FAB_eb(const Box &box,
                 const FArrayBox& src,
                 const FArrayBox &pts,
                 const int n,
                 std::string title="",
                 bool block=false);

void plot_poly_spline(const PolySpline poly, const std::string title, bool block=false);
#endif

void plot_FABs_2d(const std::map<std::string,FArrayBox>& src,
                  const int fi,
                  bool log,
                  bool block);

void plot_FABs_2d(const Vector<FArrayBox>& src,
                  const int fi,
                  std::string title="",
                  bool log=false,
                  bool block=false);

void plot_FABs_2d(const MultiFab& src,
                  const int fi,
                  std::string title="",
                  bool log=false,
                  bool block=false);

void plot_FAB_2d(const MultiFab& src,
                  const int fi,
                 const int shrink=0,
                  std::string title="",
                  bool log=false,
                  bool block=false);

void plot_FABs_2d(const iMultiFab& src,
                  const int fi,
                  std::string title="",
                  bool log=false,
                  bool block=false);



#endif // MFP_DIAGNOSTICS_H
#endif
