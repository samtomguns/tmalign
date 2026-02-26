/*
 * tmalign_core.cpp
 *
 * Stripped-down TM-align: protein-only, no file I/O, no batch processing,
 * no output formatting, no circular permutation, no mirror, no fast mode,
 * no user alignment seeding, no TMcut pre-termination.
 *
 * Exposes a single Python function:
 *
 *   tmalign(xa, ya, seqx, seqy) -> dict
 *
 * where xa and ya are (N,3) and (M,3) float64 numpy arrays of Cα coordinates,
 * seqx and seqy are single-letter amino acid strings of length N and M.
 *
 * Returns a dict with keys:
 *   tm1      - TM-score normalised by len(ya)
 *   tm2      - TM-score normalised by len(xa)
 *   rmsd     - RMSD over aligned residues within the d8 distance cutoff
 *   seq_id   - sequence identity fraction over aligned residues
 *   n_aligned - number of aligned residue pairs within d8
 *
 * Core algorithm unchanged from:
 *   Y. Zhang, J. Skolnick, Nucl. Acids Res. 33, 2302-9 (2005)
 *   Version 20220412
 */

#include <cmath>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace std;

/* -----------------------------------------------------------------------
 * Memory helpers
 * ----------------------------------------------------------------------- */

template <typename T> inline T getmin(const T &a, const T &b)
{
    return b < a ? b : a;
}

template <class A> void NewArray(A ***array, int n1, int n2)
{
    *array = new A *[n1];
    for (int i = 0; i < n1; i++) *(*array + i) = new A[n2];
}

template <class A> void DeleteArray(A ***array, int n)
{
    for (int i = 0; i < n; i++)
        if (*(*array + i)) delete[] *(*array + i);
    if (n) delete[](*array);
    (*array) = NULL;
}

/* -----------------------------------------------------------------------
 * Geometry helpers
 * ----------------------------------------------------------------------- */

// Returns squared Euclidean distance
inline double dist(double x[3], double y[3])
{
    double d0 = x[0] - y[0], d1 = x[1] - y[1], d2 = x[2] - y[2];
    return d0*d0 + d1*d1 + d2*d2;
}

inline double dot(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void transform(double t[3], double u[3][3], double *x, double *x1)
{
    x1[0] = t[0] + dot(&u[0][0], x);
    x1[1] = t[1] + dot(&u[1][0], x);
    x1[2] = t[2] + dot(&u[2][0], x);
}

void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3])
{
    for (int i = 0; i < len; i++)
        transform(t, u, &x[i][0], &x1[i][0]);
}

/* -----------------------------------------------------------------------
 * Kabsch algorithm
 * Computes optimal rotation u[3][3] and translation t[3] to superpose x→y.
 * mode=0: compute rms only; mode=1: compute t,u only; mode=2: both.
 * ----------------------------------------------------------------------- */
bool Kabsch(double **x, double **y, int n, int mode, double *rms,
            double t[3], double u[3][3])
{
    int i, j, m, m1, l, k;
    double e0, rms1, d, h, g;
    double cth, sth, sqrth, p, det, sigma;
    double xc[3], yc[3];
    double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
    const double sqrt3 = 1.73205080756888, tol = 0.01;
    const int ip[]     = {0,1,3,1,2,4,3,4,5};
    const int ip2312[] = {1,2,0,1};

    int a_failed = 0, b_failed = 0;
    const double epsilon = 1e-8;

    *rms = 0; rms1 = 0; e0 = 0;
    double c1[3], c2[3], s1[3], s2[3], sx[3], sy[3], sz[3];
    for (i = 0; i < 3; i++) { s1[i]=s2[i]=sx[i]=sy[i]=sz[i]=0.0; }
    for (i = 0; i < 3; i++) { xc[i]=yc[i]=t[i]=0.0; for (j=0;j<3;j++) { u[i][j]=0.0; r[i][j]=0.0; a[i][j]=0.0; if(i==j){u[i][j]=1.0;a[i][j]=1.0;} } }

    if (n < 1) return false;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < 3; j++) { c1[j]=x[i][j]; c2[j]=y[i][j]; s1[j]+=c1[j]; s2[j]+=c2[j]; }
        for (j = 0; j < 3; j++) { sx[j]+=c1[0]*c2[j]; sy[j]+=c1[1]*c2[j]; sz[j]+=c1[2]*c2[j]; }
    }
    for (i = 0; i < 3; i++) { xc[i]=s1[i]/n; yc[i]=s2[i]/n; }

    if (mode==2||mode==0)
        for (int mm=0;mm<n;mm++) for (int nn=0;nn<3;nn++)
            e0 += (x[mm][nn]-xc[nn])*(x[mm][nn]-xc[nn]) + (y[mm][nn]-yc[nn])*(y[mm][nn]-yc[nn]);

    for (j=0;j<3;j++) { r[j][0]=sx[j]-s1[0]*s2[j]/n; r[j][1]=sy[j]-s1[1]*s2[j]/n; r[j][2]=sz[j]-s1[2]*s2[j]/n; }

    det = r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1])
         -r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0])
         +r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);
    sigma = det;

    m = 0;
    for (j=0;j<3;j++) for (i=0;i<=j;i++) { rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j]; m++; }

    double spur = (rr[0]+rr[2]+rr[5])/3.0;
    double cof  = (((((rr[2]*rr[5]-rr[4]*rr[4])+rr[0]*rr[5])-rr[3]*rr[3])+rr[0]*rr[2])-rr[1]*rr[1])/3.0;
    det = det*det;
    for (i=0;i<3;i++) e[i]=spur;

    if (spur > 0)
    {
        d = spur*spur; h = d-cof; g = (spur*cof-det)/2.0 - spur*h;
        if (h > 0)
        {
            sqrth = sqrt(h); d = h*h*h-g*g; if(d<0.0) d=0.0;
            d = atan2(sqrt(d),-g)/3.0;
            cth = sqrth*cos(d); sth = sqrth*sqrt3*sin(d);
            e[0]=(spur+cth)+cth; e[1]=(spur-cth)+sth; e[2]=(spur-cth)-sth;

            if (mode != 0)
            {
                for (l=0; l<3; l+=2)
                {
                    d = e[l];
                    ss[0]=(d-rr[2])*(d-rr[5])-rr[4]*rr[4]; ss[1]=(d-rr[5])*rr[1]+rr[3]*rr[4];
                    ss[2]=(d-rr[0])*(d-rr[5])-rr[3]*rr[3]; ss[3]=(d-rr[2])*rr[3]+rr[1]*rr[4];
                    ss[4]=(d-rr[0])*rr[4]+rr[1]*rr[3];     ss[5]=(d-rr[0])*(d-rr[2])-rr[1]*rr[1];
                    for (int si=0;si<6;si++) if(fabs(ss[si])<=epsilon) ss[si]=0.0;

                    if (fabs(ss[0])>=fabs(ss[2])) { j=0; if(fabs(ss[0])<fabs(ss[5])) j=2; }
                    else if (fabs(ss[2])>=fabs(ss[5])) j=1;
                    else j=2;

                    d=0.0; j=3*j;
                    for (i=0;i<3;i++) { k=ip[i+j]; a[i][l]=ss[k]; d+=ss[k]*ss[k]; }
                    d = (d>epsilon) ? 1.0/sqrt(d) : 0.0;
                    for (i=0;i<3;i++) a[i][l]*=d;
                }

                d = a[0][0]*a[0][2]+a[1][0]*a[1][2]+a[2][0]*a[2][2];
                if ((e[0]-e[1])>(e[1]-e[2])) { m1=2; m=0; } else { m1=0; m=2; }
                p=0;
                for (i=0;i<3;i++) { a[i][m1]-=d*a[i][m]; p+=a[i][m1]*a[i][m1]; }
                if (p <= tol)
                {
                    p=1.0;
                    for (i=0;i<3;i++) { if(p<fabs(a[i][m])) continue; p=fabs(a[i][m]); j=i; }
                    k=ip2312[j]; l=ip2312[j+1]; p=sqrt(a[k][m]*a[k][m]+a[l][m]*a[l][m]);
                    if (p>tol) { a[j][m1]=0.0; a[k][m1]=-a[l][m]/p; a[l][m1]=a[k][m]/p; }
                    else a_failed=1;
                }
                else { p=1.0/sqrt(p); for(i=0;i<3;i++) a[i][m1]*=p; }

                if (!a_failed)
                {
                    a[0][1]=a[1][2]*a[2][0]-a[1][0]*a[2][2];
                    a[1][1]=a[2][2]*a[0][0]-a[2][0]*a[0][2];
                    a[2][1]=a[0][2]*a[1][0]-a[0][0]*a[1][2];
                }
            }
        }

        if (mode!=0 && !a_failed)
        {
            for (l=0;l<2;l++)
            {
                d=0.0;
                for (i=0;i<3;i++) { b[i][l]=r[i][0]*a[0][l]+r[i][1]*a[1][l]+r[i][2]*a[2][l]; d+=b[i][l]*b[i][l]; }
                d = (d>epsilon) ? 1.0/sqrt(d) : 0.0;
                for (i=0;i<3;i++) b[i][l]*=d;
            }
            d=b[0][0]*b[0][1]+b[1][0]*b[1][1]+b[2][0]*b[2][1]; p=0.0;
            for (i=0;i<3;i++) { b[i][1]-=d*b[i][0]; p+=b[i][1]*b[i][1]; }
            if (p<=tol)
            {
                p=1.0;
                for (i=0;i<3;i++) { if(p<fabs(b[i][0])) continue; p=fabs(b[i][0]); j=i; }
                k=ip2312[j]; l=ip2312[j+1]; p=sqrt(b[k][0]*b[k][0]+b[l][0]*b[l][0]);
                if (p>tol) { b[j][1]=0.0; b[k][1]=-b[l][0]/p; b[l][1]=b[k][0]/p; }
                else b_failed=1;
            }
            else { p=1.0/sqrt(p); for(i=0;i<3;i++) b[i][1]*=p; }

            if (!b_failed)
            {
                b[0][2]=b[1][0]*b[2][1]-b[1][1]*b[2][0];
                b[1][2]=b[2][0]*b[0][1]-b[2][1]*b[0][0];
                b[2][2]=b[0][0]*b[1][1]-b[0][1]*b[1][0];
                for (i=0;i<3;i++) for (j=0;j<3;j++)
                    u[i][j]=b[i][0]*a[j][0]+b[i][1]*a[j][1]+b[i][2]*a[j][2];
            }
            for (i=0;i<3;i++) t[i]=((yc[i]-u[i][0]*xc[0])-u[i][1]*xc[1])-u[i][2]*xc[2];
        }
    }
    else
    {
        for (i=0;i<3;i++) t[i]=((yc[i]-u[i][0]*xc[0])-u[i][1]*xc[1])-u[i][2]*xc[2];
    }

    for (i=0;i<3;i++) { if(e[i]<0) e[i]=0; e[i]=sqrt(e[i]); }
    d=e[2]; if(sigma<0.0) d=-d; d=(d+e[1])+e[0];
    if (mode==2||mode==0) { rms1=(e0-d)-d; if(rms1<0.0) rms1=0.0; }
    *rms=rms1;
    return true;
}

/* -----------------------------------------------------------------------
 * Needleman-Wunsch DP variants
 * ----------------------------------------------------------------------- */

// Variant 1: score matrix provided (used by get_initial_ssplus)
void NWDP_TM(double **score, bool **path, double **val,
             int len1, int len2, double gap_open, int j2i[])
{
    int i, j; double h, v, d;
    for (i=0;i<=len1;i++) { val[i][0]=0; path[i][0]=false; }
    for (j=0;j<=len2;j++) { val[0][j]=0; path[0][j]=false; j2i[j]=-1; }

    for (i=1;i<=len1;i++) for (j=1;j<=len2;j++)
    {
        d=val[i-1][j-1]+score[i][j];
        h=val[i-1][j]; if(path[i-1][j]) h+=gap_open;
        v=val[i][j-1]; if(path[i][j-1]) v+=gap_open;
        if (d>=h&&d>=v) { path[i][j]=true; val[i][j]=d; }
        else            { path[i][j]=false; val[i][j]=(v>=h)?v:h; }
    }
    i=len1; j=len2;
    while (i>0&&j>0)
    {
        if (path[i][j]) { j2i[j-1]=i-1; i--; j--; }
        else { h=val[i-1][j]; if(path[i-1][j]) h+=gap_open; v=val[i][j-1]; if(path[i][j-1]) v+=gap_open; if(v>=h) j--; else i--; }
    }
}

// Variant 2: distance-based with rotation (used by get_initial5 and DP_iter)
void NWDP_TM(bool **path, double **val, double **x, double **y,
             int len1, int len2, double t[3], double u[3][3],
             double d02, double gap_open, int j2i[])
{
    int i, j; double h, v, d; double xx[3], dij;
    for (i=0;i<=len1;i++) { val[i][0]=0; path[i][0]=false; }
    for (j=0;j<=len2;j++) { val[0][j]=0; path[0][j]=false; j2i[j]=-1; }

    for (i=1;i<=len1;i++)
    {
        transform(t, u, &x[i-1][0], xx);
        for (j=1;j<=len2;j++)
        {
            dij=dist(xx,&y[j-1][0]);
            d=val[i-1][j-1]+1.0/(1+dij/d02);
            h=val[i-1][j]; if(path[i-1][j]) h+=gap_open;
            v=val[i][j-1]; if(path[i][j-1]) v+=gap_open;
            if (d>=h&&d>=v) { path[i][j]=true; val[i][j]=d; }
            else            { path[i][j]=false; val[i][j]=(v>=h)?v:h; }
        }
    }
    i=len1; j=len2;
    while (i>0&&j>0)
    {
        if (path[i][j]) { j2i[j-1]=i-1; i--; j--; }
        else { h=val[i-1][j]; if(path[i-1][j]) h+=gap_open; v=val[i][j-1]; if(path[i][j-1]) v+=gap_open; if(v>=h) j--; else i--; }
    }
}

// Variant 3: secondary structure match (used by get_initial_ss)
void NWDP_TM(bool **path, double **val, const char *secx, const char *secy,
             int len1, int len2, double gap_open, int j2i[])
{
    int i, j; double h, v, d;
    for (i=0;i<=len1;i++) { val[i][0]=0; path[i][0]=false; }
    for (j=0;j<=len2;j++) { val[0][j]=0; path[0][j]=false; j2i[j]=-1; }

    for (i=1;i<=len1;i++) for (j=1;j<=len2;j++)
    {
        d=val[i-1][j-1]+1.0*(secx[i-1]==secy[j-1]);
        h=val[i-1][j]; if(path[i-1][j]) h+=gap_open;
        v=val[i][j-1]; if(path[i][j-1]) v+=gap_open;
        if (d>=h&&d>=v) { path[i][j]=true; val[i][j]=d; }
        else            { path[i][j]=false; val[i][j]=(v>=h)?v:h; }
    }
    i=len1; j=len2;
    while (i>0&&j>0)
    {
        if (path[i][j]) { j2i[j-1]=i-1; i--; j--; }
        else { h=val[i-1][j]; if(path[i-1][j]) h+=gap_open; v=val[i][j-1]; if(path[i][j-1]) v+=gap_open; if(v>=h) j--; else i--; }
    }
}

/* -----------------------------------------------------------------------
 * Parameter setters
 * ----------------------------------------------------------------------- */

void parameter_set4search(int xlen, int ylen,
    double &D0_MIN, double &Lnorm,
    double &score_d8, double &d0, double &d0_search, double &dcu0)
{
    D0_MIN = 0.5; dcu0 = 4.25;
    Lnorm  = getmin(xlen, ylen);
    d0     = (Lnorm <= 19) ? 0.168 : (1.24*pow(Lnorm-15.0, 1.0/3)-1.8);
    D0_MIN = d0 + 0.8;
    d0     = D0_MIN;
    d0_search = d0;
    if (d0_search > 8)   d0_search = 8;
    if (d0_search < 4.5) d0_search = 4.5;
    score_d8 = 1.5*pow(Lnorm, 0.3) + 3.5;
}

// Protein-only version (RNA branch removed)
void parameter_set4final(double len, double &D0_MIN, double &Lnorm,
    double &d0, double &d0_search)
{
    D0_MIN = 0.5;
    Lnorm  = len;
    d0     = (Lnorm <= 21) ? 0.5 : (1.24*pow(Lnorm-15.0, 1.0/3)-1.8);
    if (d0 < D0_MIN) d0 = D0_MIN;
    d0_search = d0;
    if (d0_search > 8)   d0_search = 8;
    if (d0_search < 4.5) d0_search = 4.5;
}

void parameter_set4scale(int len, double d_s,
    double &Lnorm, double &d0, double &d0_search)
{
    d0 = d_s; Lnorm = len; d0_search = d0;
    if (d0_search > 8)   d0_search = 8;
    if (d0_search < 4.5) d0_search = 4.5;
}

/* -----------------------------------------------------------------------
 * Score functions
 * ----------------------------------------------------------------------- */

// Normalised by Lnorm (used in search)
int score_fun8(double **xa, double **ya, int n_ali, double d,
               int i_ali[], double *score1, int score_sum_method,
               double Lnorm, double score_d8, double d0)
{
    double score_sum=0, di;
    double d_tmp=d*d, d02=d0*d0, score_d8_cut=score_d8*score_d8;
    int i, n_cut, inc=0;
    while (1)
    {
        n_cut=0; score_sum=0;
        for (i=0;i<n_ali;i++)
        {
            di=dist(xa[i],ya[i]);
            if (di<d_tmp) { i_ali[n_cut++]=i; }
            if (score_sum_method==8) { if(di<=score_d8_cut) score_sum+=1/(1+di/d02); }
            else score_sum+=1/(1+di/d02);
        }
        if (n_cut<3&&n_ali>3) { inc++; double dinc=d+inc*0.5; d_tmp=dinc*dinc; }
        else break;
    }
    *score1=score_sum/Lnorm;
    return n_cut;
}

// Normalised by n_ali (used in standard search)
int score_fun8_standard(double **xa, double **ya, int n_ali, double d,
                        int i_ali[], double *score1, int score_sum_method,
                        double score_d8, double d0)
{
    double score_sum=0, di;
    double d_tmp=d*d, d02=d0*d0, score_d8_cut=score_d8*score_d8;
    int i, n_cut, inc=0;
    while (1)
    {
        n_cut=0; score_sum=0;
        for (i=0;i<n_ali;i++)
        {
            di=dist(xa[i],ya[i]);
            if (di<d_tmp) { i_ali[n_cut++]=i; }
            if (score_sum_method==8) { if(di<=score_d8_cut) score_sum+=1/(1+di/d02); }
            else score_sum+=1/(1+di/d02);
        }
        if (n_cut<3&&n_ali>3) { inc++; double dinc=d+inc*0.5; d_tmp=dinc*dinc; }
        else break;
    }
    *score1=score_sum/n_ali;
    return n_cut;
}

/* -----------------------------------------------------------------------
 * TM-score search engine
 * Fragment-based iterative search for the best rotation among aligned pairs.
 * ----------------------------------------------------------------------- */

double TMscore8_search(double **r1, double **r2,
    double **xtm, double **ytm, double **xt, int Lali,
    double t0[3], double u0[3][3], int simplify_step, int score_sum_method,
    double *Rcomm, double local_d0_search, double Lnorm,
    double score_d8, double d0)
{
    int i, m; double score_max=-1, score, rmsd;
    const int kmax=Lali; int k_ali[kmax], ka, k;
    double t[3], u[3][3], d;
    const int n_it=20, n_init_max=6; int L_ini[n_init_max]; int L_ini_min=4;
    if (Lali<L_ini_min) L_ini_min=Lali;
    int n_init=0, i_init;
    for (i=0;i<n_init_max-1;i++) { n_init++; L_ini[i]=(int)(Lali/pow(2.0,i)); if(L_ini[i]<=L_ini_min){L_ini[i]=L_ini_min;break;} }
    if (i==n_init_max-1) { n_init++; L_ini[i]=L_ini_min; }

    int i_ali[kmax], n_cut; int L_frag, iL_max;
    for (i_init=0;i_init<n_init;i_init++)
    {
        L_frag=L_ini[i_init]; iL_max=Lali-L_frag; i=0;
        while (1)
        {
            ka=0;
            for (k=0;k<L_frag;k++) { int kk=k+i; r1[k][0]=xtm[kk][0];r1[k][1]=xtm[kk][1];r1[k][2]=xtm[kk][2]; r2[k][0]=ytm[kk][0];r2[k][1]=ytm[kk][1];r2[k][2]=ytm[kk][2]; k_ali[ka++]=kk; }
            Kabsch(r1,r2,L_frag,1,&rmsd,t,u);
            if (simplify_step!=1) *Rcomm=0;
            do_rotation(xtm,xt,Lali,t,u);
            d=local_d0_search-1;
            n_cut=score_fun8(xt,ytm,Lali,d,i_ali,&score,score_sum_method,Lnorm,score_d8,d0);
            if (score>score_max) { score_max=score; for(k=0;k<3;k++){t0[k]=t[k];u0[k][0]=u[k][0];u0[k][1]=u[k][1];u0[k][2]=u[k][2];} }

            d=local_d0_search+1;
            for (int it=0;it<n_it;it++)
            {
                ka=0;
                for (k=0;k<n_cut;k++) { m=i_ali[k]; r1[k][0]=xtm[m][0];r1[k][1]=xtm[m][1];r1[k][2]=xtm[m][2]; r2[k][0]=ytm[m][0];r2[k][1]=ytm[m][1];r2[k][2]=ytm[m][2]; k_ali[ka++]=m; }
                Kabsch(r1,r2,n_cut,1,&rmsd,t,u);
                do_rotation(xtm,xt,Lali,t,u);
                n_cut=score_fun8(xt,ytm,Lali,d,i_ali,&score,score_sum_method,Lnorm,score_d8,d0);
                if (score>score_max) { score_max=score; for(k=0;k<3;k++){t0[k]=t[k];u0[k][0]=u[k][0];u0[k][1]=u[k][1];u0[k][2]=u[k][2];} }
                if (n_cut==ka) { bool conv=true; for(k=0;k<n_cut;k++) if(i_ali[k]!=k_ali[k]){conv=false;break;} if(conv) break; }
            }
            if (i<iL_max) { i+=simplify_step; if(i>iL_max) i=iL_max; } else break;
        }
    }
    return score_max;
}

double TMscore8_search_standard(double **r1, double **r2,
    double **xtm, double **ytm, double **xt, int Lali,
    double t0[3], double u0[3][3], int simplify_step, int score_sum_method,
    double *Rcomm, double local_d0_search, double score_d8, double d0)
{
    int i, m; double score_max=-1, score, rmsd;
    const int kmax=Lali; int k_ali[kmax], ka, k;
    double t[3], u[3][3], d;
    const int n_it=20, n_init_max=6; int L_ini[n_init_max]; int L_ini_min=4;
    if (Lali<L_ini_min) L_ini_min=Lali;
    int n_init=0, i_init;
    for (i=0;i<n_init_max-1;i++) { n_init++; L_ini[i]=(int)(Lali/pow(2.0,i)); if(L_ini[i]<=L_ini_min){L_ini[i]=L_ini_min;break;} }
    if (i==n_init_max-1) { n_init++; L_ini[i]=L_ini_min; }

    int i_ali[kmax], n_cut; int L_frag, iL_max;
    for (i_init=0;i_init<n_init;i_init++)
    {
        L_frag=L_ini[i_init]; iL_max=Lali-L_frag; i=0;
        while (1)
        {
            ka=0;
            for (k=0;k<L_frag;k++) { int kk=k+i; r1[k][0]=xtm[kk][0];r1[k][1]=xtm[kk][1];r1[k][2]=xtm[kk][2]; r2[k][0]=ytm[kk][0];r2[k][1]=ytm[kk][1];r2[k][2]=ytm[kk][2]; k_ali[ka++]=kk; }
            Kabsch(r1,r2,L_frag,1,&rmsd,t,u);
            if (simplify_step!=1) *Rcomm=0;
            do_rotation(xtm,xt,Lali,t,u);
            d=local_d0_search-1;
            n_cut=score_fun8_standard(xt,ytm,Lali,d,i_ali,&score,score_sum_method,score_d8,d0);
            if (score>score_max) { score_max=score; for(k=0;k<3;k++){t0[k]=t[k];u0[k][0]=u[k][0];u0[k][1]=u[k][1];u0[k][2]=u[k][2];} }

            d=local_d0_search+1;
            for (int it=0;it<n_it;it++)
            {
                ka=0;
                for (k=0;k<n_cut;k++) { m=i_ali[k]; r1[k][0]=xtm[m][0];r1[k][1]=xtm[m][1];r1[k][2]=xtm[m][2]; r2[k][0]=ytm[m][0];r2[k][1]=ytm[m][1];r2[k][2]=ytm[m][2]; k_ali[ka++]=m; }
                Kabsch(r1,r2,n_cut,1,&rmsd,t,u);
                do_rotation(xtm,xt,Lali,t,u);
                n_cut=score_fun8_standard(xt,ytm,Lali,d,i_ali,&score,score_sum_method,score_d8,d0);
                if (score>score_max) { score_max=score; for(k=0;k<3;k++){t0[k]=t[k];u0[k][0]=u[k][0];u0[k][1]=u[k][1];u0[k][2]=u[k][2];} }
                if (n_cut==ka) { bool conv=true; for(k=0;k<n_cut;k++) if(i_ali[k]!=k_ali[k]){conv=false;break;} if(conv) break; }
            }
            if (i<iL_max) { i+=simplify_step; if(i>iL_max) i=iL_max; } else break;
        }
    }
    return score_max;
}

/* -----------------------------------------------------------------------
 * Detailed search wrappers
 * ----------------------------------------------------------------------- */

double detailed_search(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, double **x, double **y, int xlen, int ylen,
    int invmap0[], double t[3], double u[3][3], int simplify_step,
    int score_sum_method, double local_d0_search,
    double Lnorm, double score_d8, double d0)
{
    int i, j, k=0;
    for (i=0;i<ylen;i++) { j=invmap0[i]; if(j>=0){xtm[k][0]=x[j][0];xtm[k][1]=x[j][1];xtm[k][2]=x[j][2];ytm[k][0]=y[i][0];ytm[k][1]=y[i][1];ytm[k][2]=y[i][2];k++;} }
    double rmsd;
    return TMscore8_search(r1,r2,xtm,ytm,xt,k,t,u,simplify_step,score_sum_method,&rmsd,local_d0_search,Lnorm,score_d8,d0);
}

double detailed_search_standard(double **r1, double **r2,
    double **xtm, double **ytm, double **xt, double **x, double **y,
    int xlen, int ylen, int invmap0[], double t[3], double u[3][3],
    int simplify_step, int score_sum_method, double local_d0_search,
    double Lnorm, double score_d8, double d0)
{
    int i, j, k=0;
    for (i=0;i<ylen;i++) { j=invmap0[i]; if(j>=0){xtm[k][0]=x[j][0];xtm[k][1]=x[j][1];xtm[k][2]=x[j][2];ytm[k][0]=y[i][0];ytm[k][1]=y[i][1];ytm[k][2]=y[i][2];k++;} }
    double rmsd;
    return TMscore8_search_standard(r1,r2,xtm,ytm,xt,k,t,u,simplify_step,score_sum_method,&rmsd,local_d0_search,score_d8,d0);
}

/* -----------------------------------------------------------------------
 * Fast score for initialisation (3 iterations of Kabsch)
 * ----------------------------------------------------------------------- */
double get_score_fast(double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, int invmap[],
    double d0, double d0_search, double t[3], double u[3][3])
{
    double rms, tmscore, tmscore1, tmscore2; int i, j, k=0;
    for (j=0;j<ylen;j++) { i=invmap[j]; if(i>=0){ r1[k][0]=x[i][0];r1[k][1]=x[i][1];r1[k][2]=x[i][2]; r2[k][0]=y[j][0];r2[k][1]=y[j][1];r2[k][2]=y[j][2]; xtm[k][0]=x[i][0];xtm[k][1]=x[i][1];xtm[k][2]=x[i][2]; ytm[k][0]=y[j][0];ytm[k][1]=y[j][1];ytm[k][2]=y[j][2]; k++; } }
    Kabsch(r1,r2,k,1,&rms,t,u);

    int n_ali=k; double xrot[3], di;
    double d00=d0_search, d002=d00*d00, d02=d0*d0;
    double dis[n_ali]; tmscore=0;
    for (k=0;k<n_ali;k++) { transform(t,u,&xtm[k][0],xrot); di=dist(xrot,&ytm[k][0]); dis[k]=di; tmscore+=1/(1+di/d02); }

    double d002t=d002;
    while (1) { j=0; for(k=0;k<n_ali;k++) if(dis[k]<=d002t){r1[j][0]=xtm[k][0];r1[j][1]=xtm[k][1];r1[j][2]=xtm[k][2];r2[j][0]=ytm[k][0];r2[j][1]=ytm[k][1];r2[j][2]=ytm[k][2];j++;} if(j<3&&n_ali>3) d002t+=0.5; else break; }

    if (n_ali!=j) {
        Kabsch(r1,r2,j,1,&rms,t,u); tmscore1=0;
        for(k=0;k<n_ali;k++){transform(t,u,&xtm[k][0],xrot);di=dist(xrot,&ytm[k][0]);dis[k]=di;tmscore1+=1/(1+di/d02);}
        d002t=d002+1;
        while(1){j=0;for(k=0;k<n_ali;k++) if(dis[k]<=d002t){r1[j][0]=xtm[k][0];r1[j][1]=xtm[k][1];r1[j][2]=xtm[k][2];r2[j][0]=ytm[k][0];r2[j][1]=ytm[k][1];r2[j][2]=ytm[k][2];j++;} if(j<3&&n_ali>3) d002t+=0.5; else break;}
        Kabsch(r1,r2,j,1,&rms,t,u); tmscore2=0;
        for(k=0;k<n_ali;k++){transform(t,u,&xtm[k][0],xrot);di=dist(xrot,&ytm[k][0]);tmscore2+=1/(1+di/d02);}
    } else { tmscore1=tmscore; tmscore2=tmscore; }

    if (tmscore1>=tmscore) tmscore=tmscore1;
    if (tmscore2>=tmscore) tmscore=tmscore2;
    return tmscore;
}

/* -----------------------------------------------------------------------
 * Initialisation strategies
 * ----------------------------------------------------------------------- */

// Gapless threading
double get_initial(double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, int *y2x,
    double d0, double d0_search, double t[3], double u[3][3])
{
    int min_len=getmin(xlen,ylen);
    if (min_len<3) throw std::runtime_error("Sequence too short (<3)");
    int min_ali=getmin(min_len/2,5); if(min_ali<5) min_ali=5; // max(min_len/2, 5)
    // correct: min_ali = min_len/2, but at least 5
    min_ali = min_len/2; if(min_ali<=5) min_ali=5;
    int n1=-ylen+min_ali, n2=xlen-min_ali;
    int i, j, k_best=n1; double tmscore, tmscore_max=-1;
    for (int k=n1;k<=n2;k++)
    {
        for (j=0;j<ylen;j++) { i=j+k; y2x[j]=(i>=0&&i<xlen)?i:-1; }
        tmscore=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,y2x,d0,d0_search,t,u);
        if (tmscore>=tmscore_max) { tmscore_max=tmscore; k_best=k; }
    }
    for (j=0;j<ylen;j++) { i=j+k_best; y2x[j]=(i>=0&&i<xlen)?i:-1; }
    return tmscore_max;
}

// Secondary structure assignment (Cα distances only)
char sec_str(double d13, double d14, double d15, double d24, double d25, double d35)
{
    double delta=2.1;
    if (fabs(d15-6.37)<delta&&fabs(d14-5.18)<delta&&fabs(d25-5.18)<delta&&fabs(d13-5.45)<delta&&fabs(d24-5.45)<delta&&fabs(d35-5.45)<delta) return 'H';
    delta=1.42;
    if (fabs(d15-13)<delta&&fabs(d14-10.4)<delta&&fabs(d25-10.4)<delta&&fabs(d13-6.1)<delta&&fabs(d24-6.1)<delta&&fabs(d35-6.1)<delta) return 'E';
    if (d15<8) return 'T';
    return 'C';
}

void make_sec(double **x, int len, char *sec)
{
    for (int i=0;i<len;i++)
    {
        sec[i]='C';
        int j1=i-2,j2=i-1,j3=i,j4=i+1,j5=i+2;
        if (j1>=0&&j5<len)
        {
            double d13=sqrt(dist(x[j1],x[j3])), d14=sqrt(dist(x[j1],x[j4])),
                   d15=sqrt(dist(x[j1],x[j5])), d24=sqrt(dist(x[j2],x[j4])),
                   d25=sqrt(dist(x[j2],x[j5])), d35=sqrt(dist(x[j3],x[j5]));
            sec[i]=sec_str(d13,d14,d15,d24,d25,d35);
        }
    }
    sec[len]=0;
}

// SS-based alignment
void get_initial_ss(bool **path, double **val,
    const char *secx, const char *secy, int xlen, int ylen, int *y2x)
{
    NWDP_TM(path, val, secx, secy, xlen, ylen, -1.0, y2x);
}

// Local superposition
bool get_initial5(double **r1, double **r2, double **xtm, double **ytm,
    bool **path, double **val, double **x, double **y,
    int xlen, int ylen, int *y2x, double d0, double d0_search, double D0_MIN)
{
    double GL, rmsd, t[3], u[3][3];
    double d01=d0+1.5; if(d01<D0_MIN) d01=D0_MIN; double d02=d01*d01;
    double GLmax=0; int aL=getmin(xlen,ylen);
    int *invmap=new int[ylen+1];

    int n_jump1=(xlen>250)?45:(xlen>200)?35:(xlen>150)?25:15; if(n_jump1>xlen/3) n_jump1=xlen/3;
    int n_jump2=(ylen>250)?45:(ylen>200)?35:(ylen>150)?25:15; if(n_jump2>ylen/3) n_jump2=ylen/3;
    int n_frag[2]={20,100}; if(n_frag[0]>aL/3) n_frag[0]=aL/3; if(n_frag[1]>aL/2) n_frag[1]=aL/2;

    bool flag=false;
    for (int i_frag=0;i_frag<2;i_frag++)
    {
        int m1=xlen-n_frag[i_frag]+1, m2=ylen-n_frag[i_frag]+1;
        for (int i=0;i<m1;i+=n_jump1) for (int j=0;j<m2;j+=n_jump2)
        {
            for (int k=0;k<n_frag[i_frag];k++) { r1[k][0]=x[k+i][0];r1[k][1]=x[k+i][1];r1[k][2]=x[k+i][2]; r2[k][0]=y[k+j][0];r2[k][1]=y[k+j][1];r2[k][2]=y[k+j][2]; }
            Kabsch(r1,r2,n_frag[i_frag],1,&rmsd,t,u);
            NWDP_TM(path,val,x,y,xlen,ylen,t,u,d02,0.0,invmap);
            GL=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,invmap,d0,d0_search,t,u);
            if (GL>GLmax) { GLmax=GL; for(int ii=0;ii<ylen;ii++) y2x[ii]=invmap[ii]; flag=true; }
        }
    }
    delete[] invmap;
    return flag;
}

// Build score matrix combining RMSD and SS
void score_matrix_rmsd_sec(double **r1, double **r2, double **score,
    const char *secx, const char *secy, double **x, double **y,
    int xlen, int ylen, int *y2x, double D0_MIN, double d0)
{
    double t[3], u[3][3], rmsd, dij;
    double d01=d0+1.5; if(d01<D0_MIN) d01=D0_MIN; double d02=d01*d01;
    double xx[3]; int i, k=0;
    for (int j=0;j<ylen;j++) { i=y2x[j]; if(i>=0){r1[k][0]=x[i][0];r1[k][1]=x[i][1];r1[k][2]=x[i][2];r2[k][0]=y[j][0];r2[k][1]=y[j][1];r2[k][2]=y[j][2];k++;} }
    Kabsch(r1,r2,k,1,&rmsd,t,u);
    for (int ii=0;ii<xlen;ii++) {
        transform(t,u,&x[ii][0],xx);
        for (int jj=0;jj<ylen;jj++) { dij=dist(xx,&y[jj][0]); score[ii+1][jj+1]=(secx[ii]==secy[jj]) ? 1.0/(1+dij/d02)+0.5 : 1.0/(1+dij/d02); }
    }
}

// SS + local superposition combined
void get_initial_ssplus(double **r1, double **r2, double **score,
    bool **path, double **val, const char *secx, const char *secy,
    double **x, double **y, int xlen, int ylen,
    int *y2x0, int *y2x, double D0_MIN, double d0)
{
    score_matrix_rmsd_sec(r1,r2,score,secx,secy,x,y,xlen,ylen,y2x0,D0_MIN,d0);
    NWDP_TM(score,path,val,xlen,ylen,-1.0,y2x);
}

// Find longest continuous fragment (sequential Cα distances < dcu0)
void find_max_frag(double **x, int len, int *start_max, int *end_max, double dcu0)
{
    int fra_min=4, r_min=(int)(len/3.0); if(r_min>fra_min) r_min=fra_min;
    int Lfr_max=0, start=0, inc=0;
    double dcu0_cut=dcu0*dcu0, dcu_cut=dcu0_cut;
    while (Lfr_max<r_min)
    {
        Lfr_max=0; int j=1; start=0;
        for (int i=1;i<len;i++)
        {
            if (dist(x[i-1],x[i])<dcu_cut) { j++; if(i==len-1&&j>Lfr_max){Lfr_max=j;*start_max=start;*end_max=i;} }
            else { if(j>Lfr_max){Lfr_max=j;*start_max=start;*end_max=i-1;} j=1; start=i; }
        }
        if (Lfr_max<r_min) { inc++; double dinc=pow(1.1,inc)*dcu0; dcu_cut=dinc*dinc; }
    }
}

// Fragment gapless threading
double get_initial_fgt(double **r1, double **r2, double **xtm, double **ytm,
    double **x, double **y, int xlen, int ylen, int *y2x,
    double d0, double d0_search, double dcu0, double t[3], double u[3][3])
{
    int fra_min=4, fra_min1=fra_min-1;
    int xstart=0,ystart=0,xend=0,yend=0;
    find_max_frag(x,xlen,&xstart,&xend,dcu0);
    find_max_frag(y,ylen,&ystart,&yend,dcu0);
    int Lx=xend-xstart+1, Ly=yend-ystart+1;
    int L_fr=getmin(Lx,Ly);
    int *ifr=new int[L_fr], *y2x_=new int[ylen+1];

    double tmscore, tmscore_max=-1;

    // Symmetric case: both fragments same length
    if (Lx==Ly&&xlen==ylen)
    {
        int L0=xlen;
        // Part 1: use x-fragment
        for (int i=0;i<L_fr;i++) ifr[i]=xstart+i;
        if (L_fr==L0) { int n1=(int)(L0*0.1),n2=(int)(L0*0.89),j=0; for(int i=n1;i<=n2;i++) ifr[j++]=ifr[i]; L_fr=j; }
        int L1=L_fr, min_len=getmin(L1,ylen), min_ali=(int)(min_len/2.5); if(min_ali<=fra_min1) min_ali=fra_min1;
        int n1=-ylen+min_ali, n2=L1-min_ali;
        for (int k=n1;k<=n2;k++) { for(int j=0;j<ylen;j++){int i=j+k; y2x_[j]=(i>=0&&i<L1)?ifr[i]:-1;} tmscore=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,y2x_,d0,d0_search,t,u); if(tmscore>=tmscore_max){tmscore_max=tmscore;for(int j=0;j<ylen;j++) y2x[j]=y2x_[j];} }
        // Part 2: use y-fragment
        L_fr=Ly; for(int i=0;i<L_fr;i++) ifr[i]=ystart+i;
        if(L_fr==L0){int n1=(int)(L0*0.1),n2=(int)(L0*0.89),j=0;for(int i=n1;i<=n2;i++) ifr[j++]=ifr[i];L_fr=j;}
        int L2=L_fr; min_len=getmin(xlen,L2); min_ali=(int)(min_len/2.5); if(min_ali<=fra_min1) min_ali=fra_min1;
        n1=-L2+min_ali; n2=xlen-min_ali;
        for(int k=n1;k<=n2;k++){for(int j=0;j<ylen;j++) y2x_[j]=-1;for(int j=0;j<L2;j++){int i=j+k;if(i>=0&&i<xlen) y2x_[ifr[j]]=i;} tmscore=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,y2x_,d0,d0_search,t,u); if(tmscore>=tmscore_max){tmscore_max=tmscore;for(int j=0;j<ylen;j++) y2x[j]=y2x_[j];}}
        delete[] ifr; delete[] y2x_;
        return tmscore_max;
    }

    // Asymmetric case
    int L0=getmin(xlen,ylen);
    if (Lx<Ly||(Lx==Ly&&xlen<=ylen)) { for(int i=0;i<L_fr;i++) ifr[i]=xstart+i; }
    else                              { for(int i=0;i<L_fr;i++) ifr[i]=ystart+i; }
    if (L_fr==L0) { int n1=(int)(L0*0.1),n2=(int)(L0*0.89),j=0; for(int i=n1;i<=n2;i++) ifr[j++]=ifr[i]; L_fr=j; }

    if (Lx<Ly||(Lx==Ly&&xlen<=ylen))
    {
        int L1=L_fr, min_len=getmin(L1,ylen), min_ali=(int)(min_len/2.5); if(min_ali<=fra_min1) min_ali=fra_min1;
        int n1=-ylen+min_ali, n2=L1-min_ali;
        for(int k=n1;k<=n2;k++){for(int j=0;j<ylen;j++){int i=j+k;y2x_[j]=(i>=0&&i<L1)?ifr[i]:-1;} tmscore=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,y2x_,d0,d0_search,t,u); if(tmscore>=tmscore_max){tmscore_max=tmscore;for(int j=0;j<ylen;j++) y2x[j]=y2x_[j];}}
    }
    else
    {
        int L2=L_fr, min_len=getmin(xlen,L2), min_ali=(int)(min_len/2.5); if(min_ali<=fra_min1) min_ali=fra_min1;
        int n1=-L2+min_ali, n2=xlen-min_ali;
        for(int k=n1;k<=n2;k++){for(int j=0;j<ylen;j++) y2x_[j]=-1; for(int j=0;j<L2;j++){int i=j+k;if(i>=0&&i<xlen) y2x_[ifr[j]]=i;} tmscore=get_score_fast(r1,r2,xtm,ytm,x,y,xlen,ylen,y2x_,d0,d0_search,t,u); if(tmscore>=tmscore_max){tmscore_max=tmscore;for(int j=0;j<ylen;j++) y2x[j]=y2x_[j];}}
    }
    delete[] ifr; delete[] y2x_;
    return tmscore_max;
}

/* -----------------------------------------------------------------------
 * Iterative DP refinement
 * ----------------------------------------------------------------------- */
double DP_iter(double **r1, double **r2, double **xtm, double **ytm,
    double **xt, bool **path, double **val, double **x, double **y,
    int xlen, int ylen, double t[3], double u[3][3], int invmap0[],
    int g1, int g2, int iteration_max, double local_d0_search,
    double D0_MIN, double Lnorm, double d0, double score_d8)
{
    const double gap_open[2]={-0.6, 0};
    double rmsd; int *invmap=new int[ylen+1];
    int i, j, k; double tmscore, tmscore_max=-1, tmscore_old=0;
    const int score_sum_method=8, simplify_step=40;
    double d02=d0*d0;

    for (int g=g1;g<g2;g++) for (int iteration=0;iteration<iteration_max;iteration++)
    {
        NWDP_TM(path,val,x,y,xlen,ylen,t,u,d02,gap_open[g],invmap);
        k=0;
        for (j=0;j<ylen;j++) { i=invmap[j]; if(i>=0){xtm[k][0]=x[i][0];xtm[k][1]=x[i][1];xtm[k][2]=x[i][2];ytm[k][0]=y[j][0];ytm[k][1]=y[j][1];ytm[k][2]=y[j][2];k++;} }
        tmscore=TMscore8_search(r1,r2,xtm,ytm,xt,k,t,u,simplify_step,score_sum_method,&rmsd,local_d0_search,Lnorm,score_d8,d0);
        if (tmscore>tmscore_max) { tmscore_max=tmscore; for(i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
        if (iteration>0&&fabs(tmscore_old-tmscore)<1e-6) break;
        tmscore_old=tmscore;
    }
    delete[] invmap;
    return tmscore_max;
}

/* -----------------------------------------------------------------------
 * Core TM-align function (protein-only, no fast mode, no user alignment,
 * no TMcut, no circular permutation, no mirror)
 * ----------------------------------------------------------------------- */
struct TMalignResult {
    double tm_score; // TM-score under requested normalisation
    double rmsd;
    double seq_id;   // fraction of aligned pairs that are identical
    int    n_aligned; // pairs within score_d8 distance cutoff
};

TMalignResult tmalign_core(
    double **xa, double **ya,
    const char *seqx, const char *seqy,
    const char *secx, const char *secy,
    int xlen, int ylen, int norm_length)
{
    double D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0;
    double t[3], u[3][3];
    double **score, **xtm, **ytm, **xt, **r1, **r2;
    bool   **path;
    double **val;

    int minlen=getmin(xlen,ylen);
    NewArray(&score, xlen+1, ylen+1);
    NewArray(&path,  xlen+1, ylen+1);
    NewArray(&val,   xlen+1, ylen+1);
    NewArray(&xtm,   minlen, 3);
    NewArray(&ytm,   minlen, 3);
    NewArray(&xt,    xlen,   3);
    NewArray(&r1,    minlen, 3);
    NewArray(&r2,    minlen, 3);

    parameter_set4search(xlen, ylen, D0_MIN, Lnorm, score_d8, d0, d0_search, dcu0);
    const int simplify_step=40, score_sum_method=8;

    int *invmap0=new int[ylen+1];
    int *invmap =new int[ylen+1];
    for (int i=0;i<ylen;i++) invmap0[i]=-1;
    double TM, TMmax=-1;
    double ddcc = (Lnorm<=40) ? 0.1 : 0.4;
    double local_d0_search=d0_search;

    // --- Initialisation 1: gapless threading ---
    get_initial(r1,r2,xtm,ytm,xa,ya,xlen,ylen,invmap0,d0,d0_search,t,u);
    TM=detailed_search(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap0,t,u,simplify_step,score_sum_method,local_d0_search,Lnorm,score_d8,d0);
    if (TM>TMmax) TMmax=TM;
    TM=DP_iter(r1,r2,xtm,ytm,xt,path,val,xa,ya,xlen,ylen,t,u,invmap,0,2,30,local_d0_search,D0_MIN,Lnorm,d0,score_d8);
    if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }

    // --- Initialisation 2: secondary structure ---
    get_initial_ss(path,val,secx,secy,xlen,ylen,invmap);
    TM=detailed_search(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap,t,u,simplify_step,score_sum_method,local_d0_search,Lnorm,score_d8,d0);
    if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    if (TM>TMmax*0.2) {
        TM=DP_iter(r1,r2,xtm,ytm,xt,path,val,xa,ya,xlen,ylen,t,u,invmap,0,2,30,local_d0_search,D0_MIN,Lnorm,d0,score_d8);
        if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    }

    // --- Initialisation 3: local superposition ---
    if (get_initial5(r1,r2,xtm,ytm,path,val,xa,ya,xlen,ylen,invmap,d0,d0_search,D0_MIN)) {
        TM=detailed_search(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap,t,u,simplify_step,score_sum_method,local_d0_search,Lnorm,score_d8,d0);
        if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
        if (TM>TMmax*ddcc) {
            TM=DP_iter(r1,r2,xtm,ytm,xt,path,val,xa,ya,xlen,ylen,t,u,invmap,0,2,2,local_d0_search,D0_MIN,Lnorm,d0,score_d8);
            if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
        }
    }

    // --- Initialisation 4: SS + local superposition ---
    get_initial_ssplus(r1,r2,score,path,val,secx,secy,xa,ya,xlen,ylen,invmap0,invmap,D0_MIN,d0);
    TM=detailed_search(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap,t,u,simplify_step,score_sum_method,local_d0_search,Lnorm,score_d8,d0);
    if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    if (TM>TMmax*ddcc) {
        TM=DP_iter(r1,r2,xtm,ytm,xt,path,val,xa,ya,xlen,ylen,t,u,invmap,0,2,30,local_d0_search,D0_MIN,Lnorm,d0,score_d8);
        if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    }

    // --- Initialisation 5: fragment gapless threading ---
    get_initial_fgt(r1,r2,xtm,ytm,xa,ya,xlen,ylen,invmap,d0,d0_search,dcu0,t,u);
    TM=detailed_search(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap,t,u,simplify_step,score_sum_method,local_d0_search,Lnorm,score_d8,d0);
    if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    if (TM>TMmax*ddcc) {
        TM=DP_iter(r1,r2,xtm,ytm,xt,path,val,xa,ya,xlen,ylen,t,u,invmap,1,2,2,local_d0_search,D0_MIN,Lnorm,d0,score_d8);
        if (TM>TMmax) { TMmax=TM; for(int i=0;i<ylen;i++) invmap0[i]=invmap[i]; }
    }

    // Check that some alignment was found
    bool flag=false;
    for (int i=0;i<ylen;i++) if(invmap0[i]>=0){flag=true;break;}
    if (!flag) {
        delete[] invmap0; delete[] invmap;
        DeleteArray(&score,xlen+1); DeleteArray(&path,xlen+1); DeleteArray(&val,xlen+1);
        DeleteArray(&xtm,minlen); DeleteArray(&ytm,minlen); DeleteArray(&xt,xlen);
        DeleteArray(&r1,minlen); DeleteArray(&r2,minlen);
        throw std::runtime_error("No alignment found between the two structures");
    }

    // --- Final detailed search ---
    double t0[3]={0,0,0}, u0[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    TM=detailed_search_standard(r1,r2,xtm,ytm,xt,xa,ya,xlen,ylen,invmap0,t,u,1,8,local_d0_search,Lnorm,score_d8,d0);

    // Collect aligned pairs within score_d8
    int n_ali=0, k=0;
    int *m1=new int[xlen], *m2=new int[ylen];
    double d;
    do_rotation(xa,xt,xlen,t,u);
    for (int j=0;j<ylen;j++) {
        int i=invmap0[j];
        if (i>=0) {
            n_ali++;
            d=sqrt(dist(&xt[i][0],&ya[j][0]));
            if (d<=score_d8) {
                m1[k]=i; m2[k]=j;
                xtm[k][0]=xa[i][0];xtm[k][1]=xa[i][1];xtm[k][2]=xa[i][2];
                ytm[k][0]=ya[j][0];ytm[k][1]=ya[j][1];ytm[k][2]=ya[j][2];
                r1[k][0]=xt[i][0];r1[k][1]=xt[i][1];r1[k][2]=xt[i][2];
                r2[k][0]=ya[j][0];r2[k][1]=ya[j][1];r2[k][2]=ya[j][2];
                k++;
            }
        }
    }
    int n_ali8=k;

    // RMSD over aligned pairs
    double rmsd0=0;
    Kabsch(r1,r2,n_ali8,0,&rmsd0,t,u);
    rmsd0=sqrt(rmsd0/n_ali8);

    double D0_MIN2, Lnorm2, d02, d0_search2, rmsd_tmp;
    parameter_set4final((double)norm_length, D0_MIN2, Lnorm2, d02, d0_search2);
    double TM_final=TMscore8_search(r1,r2,xtm,ytm,xt,n_ali8,t0,u0,1,0,&rmsd_tmp,d0_search2,Lnorm2,score_d8,d02);

    // Sequence identity
    double Liden=0;
    for (int kk=0;kk<n_ali8;kk++) Liden+=(seqx[m1[kk]]==seqy[m2[kk]]);

    // Cleanup
    delete[] m1; delete[] m2; delete[] invmap0; delete[] invmap;
    DeleteArray(&score,xlen+1); DeleteArray(&path,xlen+1); DeleteArray(&val,xlen+1);
    DeleteArray(&xtm,minlen); DeleteArray(&ytm,minlen); DeleteArray(&xt,xlen);
    DeleteArray(&r1,minlen); DeleteArray(&r2,minlen);

    TMalignResult res;
    res.tm_score  = TM_final;
    res.rmsd      = rmsd0;
    res.seq_id    = (n_ali8>0) ? Liden/n_ali8 : 0.0;
    res.n_aligned = n_ali8;
    return res;
}

/* -----------------------------------------------------------------------
 * Python-facing wrapper
 * ----------------------------------------------------------------------- */
py::dict py_tmalign(
    py::array_t<double, py::array::c_style | py::array::forcecast> xa_in,
    py::array_t<double, py::array::c_style | py::array::forcecast> ya_in,
    const std::string &seqx,
    const std::string &seqy,
    int norm_length)
{
    auto xa_buf = xa_in.request();
    auto ya_buf = ya_in.request();

    if (xa_buf.ndim != 2 || xa_buf.shape[1] != 3)
        throw std::runtime_error("xa must be shape (N, 3)");
    if (ya_buf.ndim != 2 || ya_buf.shape[1] != 3)
        throw std::runtime_error("ya must be shape (M, 3)");

    int xlen = (int)xa_buf.shape[0];
    int ylen = (int)ya_buf.shape[0];

    if ((int)seqx.size() != xlen)
        throw std::runtime_error("seqx length must match xa.shape[0]");
    if ((int)seqy.size() != ylen)
        throw std::runtime_error("seqy length must match ya.shape[0]");
    if (xlen < 3 || ylen < 3)
        throw std::runtime_error("Both chains must have at least 3 residues");

    double *xdata = (double*)xa_buf.ptr;
    double *ydata = (double*)ya_buf.ptr;

    // Build double** view over the contiguous numpy buffers
    double **xa = new double*[xlen];
    double **ya = new double*[ylen];
    for (int i=0;i<xlen;i++) xa[i] = xdata + i*3;
    for (int i=0;i<ylen;i++) ya[i] = ydata + i*3;

    // Compute secondary structure
    char *secx = new char[xlen+1];
    char *secy = new char[ylen+1];
    make_sec(xa, xlen, secx);
    make_sec(ya, ylen, secy);

    if (norm_length < 1)
        throw std::runtime_error("norm_length must be >= 1");

    TMalignResult res = tmalign_core(xa, ya, seqx.c_str(), seqy.c_str(),
                                     secx, secy, xlen, ylen, norm_length);

    delete[] xa; delete[] ya;
    delete[] secx; delete[] secy;

    py::dict out;
    out["tm_score"]  = res.tm_score;
    out["rmsd"]      = res.rmsd;
    out["seq_id"]    = res.seq_id;
    out["n_aligned"] = res.n_aligned;
    return out;
}

/* -----------------------------------------------------------------------
 * Module definition
 * ----------------------------------------------------------------------- */
PYBIND11_MODULE(_tmalign, m)
{
    m.doc() = "TM-align: sequence-independent protein structure alignment";
    m.def("tmalign", &py_tmalign,
        py::arg("xa"), py::arg("ya"), py::arg("seqx"), py::arg("seqy"),
        py::arg("norm_length"),
        R"(
Align two protein structures and return TM-score.

Parameters
----------
xa : np.ndarray, shape (N, 3), dtype float64
    Cα coordinates of chain 1.
ya : np.ndarray, shape (M, 3), dtype float64
    Cα coordinates of chain 2.
seqx : str, length N
    Single-letter amino acid sequence of chain 1.
seqy : str, length M
    Single-letter amino acid sequence of chain 2.
norm_length : int
    Length to normalise the TM-score by (e.g. min(N, M), max(N, M), N, or M).

Returns
-------
dict with keys:
    tm_score  : TM-score normalised by norm_length
    rmsd      : RMSD over aligned residue pairs (Angstrom)
    seq_id    : sequence identity fraction over aligned pairs
    n_aligned : number of aligned residue pairs within d8 cutoff
)");
}
