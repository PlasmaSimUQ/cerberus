#ifdef PYTHON
#include "MFP_diagnostics.H"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_MultiCutFab.H>
#endif

#include "MFP_utility.H"

namespace plt = matplotlibcpp;
using namespace amrex;

void get_grid(int N, int& nr, int& nc)
{
    int n = (int)std::ceil(std::sqrt(N));

    nr = n;
    nc = n;

    for (int a=n; a>n/2; --a) {
        for (int b=n; b>n/2; --b) {
            if ((a*b >= N) && (a*b < nr*nc)) {
                nr = a;
                nc = b;
            }
        }
    }
}

void calc_gradient(const Box& box,
                   FArrayBox& data,
                   FArrayBox& slope,
                   const int d)
{

    // make sure du is empty
    slope.setVal(0.0);

    int N = data.nComp();

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array<int,3> index = {0,0,0};
    Array4<Real> const& s4 = slope.array();
    Array4<Real> const& d4 = data.array();
    //    Real dxinv;

    // cycle over dimensions

    index[d] = 1;

    for (int n=0; n<N; ++n) {
        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                        for (int i = lo.x; i <= hi.x; ++i) {
                    s4(i,j,k,n) += d4(i, j, k, n) - d4(i+index[0], j+index[1], k+index[2], n);
                }
            }
        }
    }

    return;
}

std::string print_FAB(const FArrayBox& src, int n)
{
    const Box box = src.box();
    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real const> const& src4 = src.array();

    std::stringstream s;

    s.precision(16);

    for (int k = lo.z; k <= hi.z; ++k) {
        s << "(:,:," << k << ", " << n << ")[";
        for   (int j = lo.y; j <= hi.y; ++j) {
            s << "\n[";
            for (int i = lo.x; i <= hi.x; ++i) {
                Real d = src4(i,j,0,n);
                s << d << ", ";
            }
            s << "]";
        }
        s << "]";
    }
    return s.str();
}

void plot_FAB_1d(const FArrayBox& src, std::string title, bool block)
{
    Vector<int> N(src.nComp());
    std::iota(N.begin(), N.end(), 0);

    plot_FAB_1d(src.box(), src, N, title, block);
}

void plot_FAB_1d(const FArrayBox& src, const Vector<int> N, bool block, std::string title)
{
    plot_FAB_1d(src.box(), src, N, title, block);
}


void plot_FAB_1d(const Box& box, const FArrayBox& src, const Vector<int> N, std::string title, bool block)
{

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<Real const> const& src4 = src.array();

    Vector<double> data;

    if (block)
        plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    npr = (int)std::floor(std::sqrt(N.size()));
    npc = npr + N.size() - npr*npr;

    //        npc = std::floor(std::sqrt(N.size()));
    //        npr = npc + N.size() - npc*npc;

    for (int n=0; n<N.size(); ++n) {

        data.clear();
        for (int i = lo.x; i <= hi.x; ++i) {
            data.push_back(src4(i,0,0,N[n]));
        }

        int nc = hi.x-lo.x;

        std::vector<double> xticks;

        xticks.resize(nc+1);
        std::iota (xticks.begin(), xticks.end(), lo.x);


        plt::subplot(npr,npc,n+1);
        plt::plot(xticks,data);
        plt::xlim(lo.x, hi.x);
        plt::title(std::to_string(N[n]));

    }

    plt::suptitle(title);

    if (block)
        plt::show(block);

    return;
}


#ifdef AMREX_USE_EB
PlotData2D get_2d_data(const Box& box, const EBCellFlagFab& src)
{

    PlotData2D pd;

    const Dim3 lo = amrex::lbound(box);
    const Dim3 hi = amrex::ubound(box);

    Array4<const EBCellFlag> const& src4 = src.array();

    Vector<float> &data = pd.data;

    for   (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            float d = (float)src4(i,j,0).isCovered();
            data.push_back(d);
        }
    }

    pd.nr = hi.y-lo.y + 1;
    pd.nc = hi.x-lo.x + 1;

    pd.extents = {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f};

    return pd;
}


void plot_FAB_2d(const Box& box, const EBCellFlagFab& src, std::string title, bool block)
{

    plt::figure_size(1200, 780);

    PlotData2D pd = get_2d_data(box, src);


    plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);
    plt::colorbar();


    plt::title(title);


    if (block)
        plt::show(block);

    return;
}

void plot_FAB_2d(const EBCellFlagFab& src, std::string title, bool block)
{
    plot_FAB_2d(src.box(), src, title, block);
    return;
}

void plot_FAB_eb(const Box &box, const FArrayBox& src, const FArrayBox &pts, const int n, std::string title, bool block)
{
    if (block)
        plt::figure_size(1200, 780);

    Array4<const Real> const& pts4 = pts.array();

    // plot gradient
    PlotData2D pd = get_2d_data(src.box(), src, n);
    plt::imshow(pd.data.data(), pd.nr, pd.nc, 1, {{"origin","lower"}}, pd.extents);
    plt::colorbar();

    const Box &pts_box = pts.box();
    Dim3 lo = amrex::lbound(pts_box);
    Dim3 hi = amrex::ubound(pts_box);

    for     (int k = lo.z; k <= hi.z; ++k) {
        for   (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
                    for (int i = lo.x; i <= hi.x; ++i) {
                // plot cells

                Vector<double> box_x = {-0.5, 0.5, 0.5, -0.5, -0.5};
                Vector<double> box_y = {-0.5, -0.5, 0.5, 0.5, -0.5};

                for (int bi=0; bi<5; ++bi) {
                    box_x[bi] += i;
                    box_y[bi] += j;
                }

                std::map<std::string,std::string> string_args = {{"color","k"}, {"linestyle","dotted"}};
                std::map<std::string,double> numeric_args = {{"alpha",0.3}};
                plt::plot(box_x, box_y, string_args, numeric_args);

                int cnt = (int) pts4(i,j,k,0);

                if (cnt < 2) continue;

                // plot the poly
                Vector<double> plot_x, plot_y;


                for (int n=0; n<cnt; ++n) {
                    plot_x.push_back(pts4(i,j,k,n*AMREX_SPACEDIM + 1) + i);
                    plot_y.push_back(pts4(i,j,k,n*AMREX_SPACEDIM + 2) + j);
                }

                plt::plot(plot_x, plot_y, "k");

            }
        }
    }

    // bounding box
    lo = amrex::lbound(box);
    hi = amrex::ubound(box);

    Vector<double> box_x = {-0.5+lo.x, 0.5+hi.x, 0.5+hi.x, -0.5+lo.x, -0.5+lo.x};
    Vector<double> box_y = {-0.5+lo.y, -0.5+lo.y, 0.5+hi.y, 0.5+hi.y, -0.5+lo.y};

    plt::plot(box_x, box_y,"r--");

    plt::suptitle(title);

    if (block)
        plt::show(block);
}

void plot_FAB_eb(const FArrayBox& src, const FArrayBox &pts, const int n, std::string title, bool block) {
    plot_FAB_eb(src.box(), src, pts, n, title, block);
}

#endif

void plot_FABs_2d(const std::map<std::string,FArrayBox>& src, const int fi, bool log, bool block)
{

    if (block)
        plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    int N = src.size();

    get_grid(N, npr, npc);

//    npr = (int)std::floor(std::sqrt(N));
//    npc = npr + N - npr*npr;

    //        npc = std::floor(std::sqrt(N));
    //        npr = npc + N - npc*npc;

    int n = 1;
    for (const auto& FAB : src) {

        const Box& box = FAB.second.box();

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        Array4<Real const> const& src4 = FAB.second.array();

        Vector<float> data;

        data.clear();
        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float d = (float)src4(i,j,0,fi);
                if (log) {
                    d = std::log10(std::abs(d));
                }
                data.push_back(d);
            }
        }

        int nr = hi.y-lo.y;
        int nc = hi.x-lo.x;

        std::vector<int> xticks, yticks;
        std::vector<std::string> xlabels, ylabels;

        xticks.resize(nc+1);
        std::iota (xticks.begin(), xticks.end(), lo.x);
        for (const auto &i : xticks) {
            xlabels.push_back(std::to_string(i));
        }

        yticks.resize(nr+1);
        std::iota (yticks.begin(), yticks.end(), lo.y);
        for (const auto &i : yticks) {
            ylabels.push_back(std::to_string(i));
        }


        plt::subplot(npr,npc,n);
        plt::title(FAB.first);
        plt::imshow(data.data(), nr+1, nc+1, 1, {{"origin","lower"}}, {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f});
        //        plt::pcolormesh(xlabel, ylabel, data,{{"edgecolor", "k"}});
        plt::colorbar();

        n++;
    }

//    plt::suptitle(title);
    plt::tight_layout();

    if (block)
        plt::show(block);

    return;
}

void plot_FABs_2d(const Vector<FArrayBox>& src, const int fi, std::string title, bool log, bool block)
{

    if (block)
        plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    int N = src.size();

    npr = (int)std::floor(std::sqrt(N));
    npc = npr + N - npr*npr;

    //        npc = std::floor(std::sqrt(N));
    //        npr = npc + N - npc*npc;

    int n = 1;
    for (const auto& FAB : src) {

        const Box& box = FAB.box();

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        Array4<Real const> const& src4 = FAB.array();

        Vector<float> data;

        data.clear();
        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float d = (float)src4(i,j,0,fi);
                if (log) {
                    d = std::log10(std::abs(d));
                }
                data.push_back(d);
            }
        }

        int nr = hi.y-lo.y;
        int nc = hi.x-lo.x;

        std::vector<int> xticks, yticks;
        std::vector<std::string> xlabels, ylabels;

        xticks.resize(nc+1);
        std::iota (xticks.begin(), xticks.end(), lo.x);
        for (const auto &i : xticks) {
            xlabels.push_back(std::to_string(i));
        }

        yticks.resize(nr+1);
        std::iota (yticks.begin(), yticks.end(), lo.y);
        for (const auto &i : yticks) {
            ylabels.push_back(std::to_string(i));
        }


        plt::subplot(npr,npc,n);
        plt::imshow(data.data(), nr+1, nc+1, 1, {{"origin","lower"}}, {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f});
        //        plt::pcolormesh(xlabel, ylabel, data,{{"edgecolor", "k"}});
        plt::colorbar();

        n++;
    }

    plt::suptitle(title);

    if (block)
        plt::show(block);

    return;
}

void plot_FABs_2d(const MultiFab& src, const int fi, std::string title, bool log, bool block)
{

    if (block)
        plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    int N = src.size();

    npr = (int)std::floor(std::sqrt(N));
    npc = npr + N - npr*npr;

    //        npc = std::floor(std::sqrt(N));
    //        npr = npc + N - npc*npc;

    int n = 1;
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {

        const FArrayBox &FAB = src[mfi];

        const Box& box = FAB.box();

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        Array4<Real const> const& src4 = FAB.array();

        Vector<float> data;

        data.clear();
        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float d = (float)src4(i,j,0,fi);
                if (log) {
                    d = std::log10(std::abs(d));
                }
                data.push_back(d);
            }
        }

        int nr = hi.y-lo.y;
        int nc = hi.x-lo.x;

        std::vector<int> xticks, yticks;
        std::vector<std::string> xlabels, ylabels;

        xticks.resize(nc+1);
        std::iota (xticks.begin(), xticks.end(), lo.x);
        for (const auto &i : xticks) {
            xlabels.push_back(std::to_string(i));
        }

        yticks.resize(nr+1);
        std::iota (yticks.begin(), yticks.end(), lo.y);
        for (const auto &i : yticks) {
            ylabels.push_back(std::to_string(i));
        }


        plt::subplot(npr,npc,n);
        plt::imshow(data.data(), nr+1, nc+1, 1, {{"origin","lower"}}, {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f});
        //        plt::pcolormesh(xlabel, ylabel, data,{{"edgecolor", "k"}});
        plt::title(num2str(mfi.index()));
        plt::colorbar();

        n++;
    }

    plt::suptitle(title);

    if (block)
        plt::show(block);

    return;
}


void plot_FAB_2d(const MultiFab& src, const int fi, const int shrink, std::string title, bool log, bool block)
{

    plt::figure_size(1200, 780);



    int N = src.size();

    // figure out subplot array sizes
//    int npr, npc;
//    npr = (int)std::floor(std::sqrt(N));
//    npc = npr + N - npr*npr;

    //        npc = std::floor(std::sqrt(N));
    //        npr = npc + N - npc*npc;

    // get the size of the box we need to plot
    int xlo = std::numeric_limits<int>::max();
    int xhi = -xlo;
    int ylo = xlo;
    int yhi = xhi;

    for (MFIter mfi(src); mfi.isValid(); ++mfi) {

        const FArrayBox& FAB = src[mfi];
        const Box& box = FAB.box();

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        xlo = std::min(xlo, lo.x);
        xhi = std::max(xhi, hi.x+1);

        ylo = std::min(ylo, lo.y);
        yhi = std::max(yhi, hi.y+1);
    }

    int nc = xhi-xlo;
    int nr = yhi-ylo;

    Vector<float> data(nc*nr);

    std::vector<int> xticks, yticks;
    std::vector<std::string> xlabels, ylabels;

    xticks.resize(nc+1);
    std::iota (xticks.begin(), xticks.end(), xlo);
    for (const auto &i : xticks) {
        xlabels.push_back(std::to_string(i));
    }

    yticks.resize(nr+1);
    std::iota (yticks.begin(), yticks.end(), ylo);
    for (const auto &i : yticks) {
        ylabels.push_back(std::to_string(i));
    }

    for (MFIter mfi(src); mfi.isValid(); ++mfi) {

        const FArrayBox& FAB = src[mfi];

        const Box box = grow(FAB.box(), -shrink);


        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        xlo = std::min(xlo, lo.x);
        xhi = std::max(xhi, hi.x+1);

        ylo = std::min(ylo, lo.y);
        yhi = std::max(yhi, hi.y+1);

        Array4<Real const> const& src4 = FAB.array();

        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float d = (float)src4(i,j,0,fi);
                if (log) {
                    d = std::log10(std::abs(d));
                }
                data[(j-ylo)*nc+(i-xlo)] = d;
            }
        }
    }

    plt::imshow(data.data(), nr, nc, 1, {{"origin","lower"}}, {xlo-0.5f, xhi+0.5f, ylo-0.5f, yhi+0.5f});


    plt::colorbar();
    plt::suptitle(title);

    if (block)
        plt::show(block);

    return;
}


void plot_FABs_2d(const iMultiFab& src, const int fi, std::string title, bool log, bool block)
{

    if (block)
        plt::figure_size(1200, 780);

    // figure out subplot array sizes
    int npr, npc;

    int N = src.size();

    npr = (int)std::floor(std::sqrt(N));
    npc = npr + N - npr*npr;

    //        npc = std::floor(std::sqrt(N));
    //        npr = npc + N - npc*npc;

    int n = 1;

    for (MFIter mfi(src); mfi.isValid(); ++mfi) {

        const auto& FAB = src[mfi];

        const Box& box = FAB.box();

        const Dim3 lo = amrex::lbound(box);
        const Dim3 hi = amrex::ubound(box);

        auto const& src4 = FAB.array();

        Vector<float> data;

        data.clear();
        for   (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                float d = (float)src4(i,j,0,fi);
                if (log) {
                    d = std::log10(std::abs(d));
                }
                data.push_back(d);
            }
        }

        int nr = hi.y-lo.y;
        int nc = hi.x-lo.x;

        std::vector<int> xticks, yticks;
        std::vector<std::string> xlabels, ylabels;

        xticks.resize(nc+1);
        std::iota (xticks.begin(), xticks.end(), lo.x);
        for (const auto &i : xticks) {
            xlabels.push_back(std::to_string(i));
        }

        yticks.resize(nr+1);
        std::iota (yticks.begin(), yticks.end(), lo.y);
        for (const auto &i : yticks) {
            ylabels.push_back(std::to_string(i));
        }


        plt::subplot(npr,npc,n);
        plt::imshow(data.data(), nr+1, nc+1, 1, {{"origin","lower"}}, {lo.x-0.5f, hi.x+0.5f, lo.y-0.5f, hi.y+0.5f});
        //        plt::pcolormesh(xlabel, ylabel, data,{{"edgecolor", "k"}});
        plt::colorbar();

        n++;
    }

    plt::suptitle(title);

    if (block)
        plt::show(block);

    return;
}

#endif
