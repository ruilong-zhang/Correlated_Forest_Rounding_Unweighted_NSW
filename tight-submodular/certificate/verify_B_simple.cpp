// Standalone directed cover of B for the compressed tight NSW menu.
// Proof: menu_reductions.tex and uniform_middle_band.tex.
// B omits redundant NEW sqrt and physical private-guard constraints.
// All proof decisions use exact integer arithmetic; floating point chooses
// only Fenchel tangents and the order of bisections.  Assertions are required.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
using i64=std::int64_t;
using i128=__int128_t;
constexpr int BITS=48;
constexpr i64 SCALE=i64(1)<<BITS, GUARD=i64(1)<<62;
i64 checked(i128 x){assert(-i128(GUARD)<x && x<i128(GUARD));return i64(x);}
i64 down(i128 a,i128 b){assert(b>0);i128 q=a/b,r=a%b;if(r<0)--q;return checked(q);}
i64 up(i128 a,i128 b){assert(b>0);i128 q=a/b,r=a%b;if(r>0)++q;return checked(q);}
struct I {i64 l,h;I(i64 a=0,i64 b=0):l(a),h(b){assert(l<=h);checked(l);checked(h);}};
I pt(i64 a){return I(a,a);} I num(i64 a){return pt(checked(i128(a)*SCALE));}
I rat(i64 a,i64 b){return I(down(i128(a)*SCALE,b),up(i128(a)*SCALE,b));}
I operator+(I a,I b){return I(checked(i128(a.l)+b.l),checked(i128(a.h)+b.h));}
I operator-(I a){return I(-a.h,-a.l);} I operator-(I a,I b){return a+-b;}
I operator*(I a,I b){i128 v[]={i128(a.l)*b.l,i128(a.l)*b.h,i128(a.h)*b.l,i128(a.h)*b.h};auto z=std::minmax_element(v,v+4);return I(down(*z.first,SCALE),up(*z.second,SCALE));}
I operator/(I a,I b){assert(b.l>0);i64 l=GUARD-1,h=-GUARD+1;for(i64 x:{a.l,a.h})for(i64 y:{b.l,b.h}){l=std::min(l,down(i128(x)*SCALE,y));h=std::max(h,up(i128(x)*SCALE,y));}return I(l,h);}
I mx(I a,I b){return I(std::max(a.l,b.l),std::max(a.h,b.h));}
I mn(I a,I b){return I(std::min(a.l,b.l),std::min(a.h,b.h));}
I clip(I a,I lo,I hi){return mn(mx(a,lo),hi);}
I hull(I a,I b){return I(std::min(a.l,b.l),std::max(a.h,b.h));}
I ZERO=num(0),ONE=num(1);
long double fl(I a){return (static_cast<long double>(a.l)+a.h)/(2*SCALE);}
std::vector<I> LOG;
std::array<std::vector<I>,20> EPS;
std::map<std::string,I> CONST;
int LOG_N=0,EPS_N=0;
I L,Z,BETA,KAPPA,ZETA,RSTAR,HOME_JOIN,MUR,T0,B,KM,ETA,LAM,U,C0;
// Directed literals independently enclosed in check_localized_home_credit.py.
I EXP_A(281756592471775,281756592471783);
I CAP_B(280072281368,280072281376);
I HALF_LOG_R4(45531063257,45531063265);
I HALF_PRICE_CAP=rat(7,6)*HALF_LOG_R4;
// B pays b=min(sigma,Pc)/3 in [0,Pc/3], Pc=(43/25) log R4.
I PRIVATE_PAYMENT_CAP=rat(86,75)*HALF_LOG_R4;
I ALOCAL=rat(1,1000);
I readI(std::istream& f){i64 a,b;f>>a>>b;assert(f);return I(a,b);}
void token(std::istream& f,const std::string& want){std::string got;f>>got;assert(f && got==want);}

I logpoint(i64 x){
 assert(x>0);int e=0;i64 base=SCALE;
 while(x<base){base>>=1;--e;assert(base>0);}
 while(i128(x)>=2*i128(base)){base=checked(i128(base)*2);++e;}
 i128 p=i128(x-base)*LOG_N;i64 j=down(p,base);assert(j>=0 && j<LOG_N);
 I w=I(down((p-i128(j)*base)*SCALE,base),up((p-i128(j)*base)*SCALE,base));
 I lower=pt(LOG[j].l)+w*(pt(LOG[j+1].l)-pt(LOG[j].l));
 I m=I(down(i128(x)*SCALE,base),up(i128(x)*SCALE,base));
 I x0=ONE+rat(j,LOG_N),x1=ONE+rat(j+1,LOG_N);
 I upper=mn(pt(LOG[j].h)+(m-x0)/x0,pt(LOG[j+1].h)+(m-x1)/x1);
 return I(lower.l,upper.h)+num(e)*L;
}
I lg(I x){assert(x.l>0);return I(logpoint(x.l).l,logpoint(x.h).h);}

I eps_upper(I r){
 assert(r.l>=0);
 const i64 first=SCALE>>10,last=SCALE<<10;
 if(r.l<first)return I(0,SCALE);
 if(r.l>=last)return I(0,(ONE/(num(8)*pt(r.l))).h);
 i64 base=first;int e=-10;
 while(i128(r.l)>=2*i128(base)){base*=2;++e;}
 i128 p=i128(r.l-base)*EPS_N;i64 j=down(p,base);assert(j>=0 && j<EPS_N);
 I w=I(down((p-i128(j)*base)*SCALE,base),up((p-i128(j)*base)*SCALE,base));
 I v=pt(EPS[e+10][j].h)+w*(pt(EPS[e+10][j+1].h)-pt(EPS[e+10][j].h));
 return I(0,v.h);
}

struct Proxy{I tau,alpha,k,R;};
struct Planar{I a,c,d;};
std::vector<Proxy> PROX;
std::vector<Planar> PLANAR;
std::vector<I> PROFILE;
std::vector<std::vector<i64>> PROFILE_RMQ;
I TMAX,DMAX;
void read_tables(const std::string& path){
 std::ifstream f(path);assert(f);int n,p;token(f,"P");f>>p;assert(p==BITS);
 token(f,"LOG_N");f>>LOG_N;token(f,"LOG");f>>n;assert(n==LOG_N+1);
 for(int j=0;j<n;++j)LOG.push_back(readI(f));
 token(f,"EPS_N");f>>EPS_N;token(f,"EPS_BLOCKS");f>>n;assert(n==20);
 for(int k=0;k<20;++k){int e,nn;token(f,"EPS");f>>e>>nn;assert(e==k-10 && nn==EPS_N+1);for(int j=0;j<nn;++j)EPS[k].push_back(readI(f));}
 token(f,"CONSTANTS");f>>n;for(int j=0;j<n;++j){std::string s;f>>s;CONST[s]=readI(f);}token(f,"END");
 L=CONST.at("L");Z=CONST.at("Z");BETA=CONST.at("beta");KAPPA=CONST.at("kappa");
 ZETA=CONST.at("zeta");RSTAR=CONST.at("rstar");HOME_JOIN=CONST.at("home_join");
 B=rat(9,16);KM=rat(5,16);ETA=rat(3,5);LAM=rat(2319048039,10000000000);U=rat(6905799,10000000);C0=rat(3999999,10000000); MUR=ONE/(num(4)*BETA); T0=(ONE+lg(MUR)-PRIVATE_PAYMENT_CAP)/MUR;
 for(auto kv:std::vector<std::pair<std::string,I>>{{"b",B},{"k_M",KM},{"eta",ETA},{"c0",C0}}){I v=CONST.at(kv.first);assert(v.l<=kv.second.l && v.h>=kv.second.h);}
}
void read_input(const std::string& path){
 std::ifstream f(path);assert(f);int n,p;token(f,"P");f>>p;assert(p==BITS);
 token(f,"PROXIES");f>>n;for(int j=0;j<n;++j)PROX.push_back({readI(f),readI(f),readI(f),readI(f)});
 token(f,"PLANAR");f>>n;for(int j=0;j<n;++j)PLANAR.push_back({readI(f),readI(f),readI(f)});
 token(f,"PROFILE");f>>n;assert(n==4001);for(int j=0;j<n;++j)PROFILE.push_back(readI(f));token(f,"END");
 TMAX=ONE/C0;DMAX=ONE/PROFILE[0];assert(TMAX.h<rat(2500000626,1000000000).l && DMAX.h<num(21).l);
 // At m=j/1000 the implemented continuous curve equals G_(j-1).
 PROFILE_RMQ.push_back(std::vector<i64>(4001,GUARD-1));
 for(int j=1;j<=4000;++j)PROFILE_RMQ[0][j]=(pt(PROFILE[j-1].l)/rat(j,1000)).l;
 for(int d=1;(1<<d)<=4001;++d){int len=4001-(1<<d)+1;PROFILE_RMQ.emplace_back(len);for(int j=0;j<len;++j)PROFILE_RMQ[d][j]=std::min(PROFILE_RMQ[d-1][j],PROFILE_RMQ[d-1][j+(1<<(d-1))]);}
}

I M(I v){assert(v.l>=0);I x=B*v;return lg(mn(ONE+x*x/(ONE+x),ONE+KM*v));}
I home_logW(I T,I D){
 if(D.h==0)return lg(mx(ONE,T));
 assert(D.l>0);
 I flat=lg(mx(ONE,D));
 if(T.h<=(BETA*D).l)return flat;
 if(T.l>=(HOME_JOIN*D).h)return lg(mx(ONE,T+ZETA*D));
 I r=T/D-BETA;
 if(r.h<=0)return flat;
 I exponential=mx(ZERO,lg(D)+L*mx(ZERO,r));
 if(r.h<=RSTAR.l)return exponential;
 I affine=lg(mx(ONE,T+ZETA*D));
 if(r.l>=RSTAR.h)return affine;
 return hull(exponential,affine);
}
I home_onset(I D){
 if(D.l>=SCALE)return BETA*D;
 if(D.h<=L.l)return ONE-ZETA*D;
 if(D.l>=L.h)return D*(BETA-lg(D)/L);
 return hull(ONE-ZETA*D,D*(BETA-lg(D)/L));
}
I f10_small(I s){
 assert(s.l>=0 && s.h<=ALOCAL.h);
 // In f10(s)/s the term of order x^k, x=10s, has sign (-1)^(k+1).
 // Thus quotient order 7 is an upper polynomial, and order 8 lower.
 I x=num(10)*s,term=s*x/num(2),sum=term,upper=sum;
 for(int k=2;k<=8;++k){term=term*x/num(k+1);sum=(k%2)?sum+term:sum-term;if(k==7)upper=sum;}
 return mx(ZERO,I(sum.l,upper.h));
}
I localized_h_point(i64 s){
 assert(s>=0);I x=pt(s);
 if(s<=ALOCAL.l)return f10_small(x);
 if(s>=ALOCAL.h)return mx(ZERO,x-CAP_B);
 return mx(ZERO,hull(f10_small(x),x-CAP_B));
}
I localized_h(I s){
 s=mx(ZERO,s);return I(localized_h_point(s.l).l,localized_h_point(s.h).h);
}
I home_level_a(I D){
 assert(D.l>0);
 I linear=EXP_A-ZETA*D;
 I exponential=D*(BETA+(ALOCAL-lg(D))/L);
 if(D.h<=(L*EXP_A).l)return linear;
 if(D.l>=(L*EXP_A).h)return exponential;
 return hull(linear,exponential);
}
I home_C(I t,I D,I q,I mu){
 assert(t.l==t.h && D.l==D.h && q.l==q.h && mu.l==mu.h);
 assert(mu.l>=rat(1,100).h);
 if(D.h==0)return mu*q;
 I qt=Z*t;q=mn(q,pt(qt.l));
 I qs[5]={q,qt,clip(t-home_onset(D),q,qt),clip(t-HOME_JOIN*D,q,qt),clip(t-home_level_a(D),q,qt)};
 I md=M(D),ans=num(1000);
 bool inactive=D.l>=rat(401,400).h || (BETA*t+ZETA*D).l>=rat(401,400).h;
 int count=inactive?4:5;
 for(int j=0;j<count;++j){I Q=qs[j],T=mx(ZERO,t-Q),s=home_logW(T,D);I credit=mn(md,ETA*s);if(!inactive)credit=mn(credit,localized_h(s));ans=mn(ans,mu*Q+credit);}
 return ans;
}

// Evaluate the implemented delayed private curve at one dyadic point.
I private_profile(I t,I D){
 if(D.h==0)return ZERO;
 assert(t.l==t.h && D.l==D.h);
 if(i128(t.l)*1000>=i128(D.l)*4000)return D*PROFILE[3999];
 i64 j=down(i128(t.l)*1000,D.l);assert(j>=0 && j<4000);
 I x=t/D*num(1000)-num(j);
 I lo=PROFILE[std::max<i64>(0,j-1)],hi=PROFILE[j];
 return D*(lo+x*(hi-lo));
}
i64 range_profile_ratio(int lo,int hi){
 assert(1<=lo && lo<=hi && hi<=4000);int d=0;while((1<<(d+1))<=hi-lo+1)++d;
 return std::min(PROFILE_RMQ[d][lo],PROFILE_RMQ[d][hi-(1<<d)+1]);
}
I private_profile_lower(I t,I D){
 // Monotonicity in t reduces to t.l.  As a function of D the score is
 // piecewise affine; its extrema are at endpoints or m-grid crossings.
 I ans=mn(private_profile(pt(t.l),pt(D.l)),private_profile(pt(t.l),pt(D.h)));
 if(t.l==0 || D.h==0)return ans;
 i64 first=up(i128(t.l)*1000,D.h),last=D.l?down(i128(t.l)*1000,D.l):4000;
 int l=int(std::max<i64>(1,first)),h=int(std::min<i64>(4000,last));
 if(l<=h)ans=mn(ans,pt(t.l)*pt(range_profile_ratio(l,h)));
 return ans;
}
bool private_excluded(I t,I D){
 if((C0*pt(t.l)).l>SCALE)return true;
 if(private_profile_lower(t,D).l>SCALE)return true;
 if(t.l>D.h && t.l>SCALE){I test=num(8)*(pt(t.l)-ONE)*(pt(t.l)-pt(D.h))-pt(D.h)*pt(D.h);if(test.l>0)return true;}
 return false;
}

struct Key{std::array<i64,4> v;bool operator==(const Key& b)const{return v==b.v;}};
struct Hash{std::size_t operator()(const Key& a)const{std::size_t h=0;for(i64 x:a.v)h^=std::hash<i64>{}(x)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
struct TDInfo{I mu;bool excluded;};
std::unordered_map<Key,TDInfo,Hash> TD_CACHE;
I mu_upper(I t){return t.h<=SCALE?ONE:mn(ONE,ONE/t);}
I root_function(I tC,I DC,I t_cost,I mu){return home_C(tC,DC,ZERO,mu)-mu*t_cost+ONE+lg(mu);}
TDInfo td_info(I t,I D){
 Key key{{t.l,t.h,D.l,D.h}};auto found=TD_CACHE.find(key);if(found!=TD_CACHE.end())return found->second;
 TDInfo out;out.excluded=private_excluded(t,D);
 if(out.excluded){out.mu=MUR;TD_CACHE.emplace(key,out);return out;}
 if(t.h<=T0.l && D.h<=SCALE)out.mu=MUR;
 else{
  I upper=I(mu_upper(pt(t.h)).l,mu_upper(pt(t.l)).h);
  I floor=MUR/mx(ONE,D);
  assert(floor.l>=rat(1,100).h);
  // Search the clipped multiplier directly, never invoking the five-point
  // credit lemma at a trial below the global positive floor.
  i64 lo=floor.l,hi=upper.h;
  // A lower enclosure of the root: f_actual <= f_upper.
  for(int j=0;j<31 && hi-lo>1;++j){i64 mid=lo+(hi-lo)/2;I f=root_function(pt(t.h),pt(D.h),pt(t.l),pt(mid));if(f.h<0)lo=mid;else hi=mid;}
  i64 rootlo=lo,rootup=upper.h;
  i64 candidate=mu_upper(pt(t.h)).l;
  I ftop=root_function(pt(t.l),pt(D.l),pt(t.h),pt(candidate))-PRIVATE_PAYMENT_CAP;
  if(ftop.l>=0){
   lo=floor.l;hi=std::max(candidate,lo);
   for(int j=0;j<31 && hi-lo>1;++j){i64 mid=lo+(hi-lo)/2;I f=root_function(pt(t.l),pt(D.l),pt(t.h),pt(mid))-PRIVATE_PAYMENT_CAP;if(f.l>=0)hi=mid;else lo=mid;}
   rootup=hi;
  }
  out.mu=mn(upper,mx(I(rootlo,rootup),floor));
  assert(out.mu.l>0 && out.mu.h<=SCALE+1);
 }
 if(TD_CACHE.size()>500000)TD_CACHE.clear();
 TD_CACHE.emplace(key,out);return out;
}

I gap_upper(I G,I floor){
 assert(G.l>=0 && floor.l>=0);
 if(G.h==0)return ZERO;
 if(G.l<1024 || floor.l==0)return I(0,G.h);
 if(i128(floor.l)>=i128(1024)*G.h)return I(0,(G*G/(num(8)*floor)).h);
 return G*eps_upper(floor/G);
}
I joint_score(const Planar& p,I t,I D,I y,I ss){
 I total=p.c+p.d,mean=(p.a+p.d)*y+total*t-p.d*ss,width=total*D;
 I floor=mx(ZERO,mx((p.a+p.d)*y,mean-width));
 return mean-gap_upper(width,floor);
}
// Exactly eight terms: 1, singleton, five affine rows, one joint row.
// The five affine rows form two exact envelopes of stored functions.
// On B, OLD .4/R4 dominates every eligible NEW sqrt score; c0*t<=1
// follows from the private feasibility test, which is retained unchanged.
int score_count(){return 2+int(PROX.size()+PLANAR.size());}
I physical_score(int index,I t,I D,I y,I ss){
 if(index==0)return ONE;
 if(index==1)return y;
 index-=2;
 if(index<int(PROX.size())){const auto& p=PROX[index];return (y+p.k*t-(ONE-p.alpha)*p.k*ss)/p.R;}
 return joint_score(PLANAR[index-PROX.size()],t,D,y,ss);
}
bool score_active(int index,I D,I y){
 index-=2;if(index<0 || index>=int(PROX.size()))return true;
 return PROX[index].tau.h==0 || D.h<=(PROX[index].tau*pt(y.l)).l;
}
I physical(I t,I D,I y,I ss){
 I ans=ONE;for(int j=0;j<score_count();++j)if(score_active(j,D,y))ans=mx(ans,physical_score(j,t,D,y,ss));return ans;
}
I localized_g_point(i64 w){
 assert(w>0);I v=pt(w);
 if(w<=SCALE)return v-ONE;
 I inv=ONE/v,pow=ONE;for(int j=0;j<10;++j)pow=pow*inv;
 return mn((ONE-pow)/num(10),CAP_B);
}
I localized_g(I w){return I(localized_g_point(w.l).l,localized_g_point(w.h).h);}
I owner_net(I w,I y,I D){
 I share=(y.l==0 || D.h>=y.l)?ONE:D/y;
 I logw=lg(w),old=logw-mn(M(w*share),ETA*logw);
 return mx(ZERO,mx(old-HALF_PRICE_CAP,rat(499,1000)*localized_g(w)));
}

int H=1;
struct Box{I t,D;std::array<I,4> y,ss;int depth=0;};
bool propagate(Box& b){
 b.t.h=std::min(b.t.h,TMAX.h);b.D.h=std::min(b.D.h,DMAX.h);
 for(int j=0;j<H;++j){b.t.l=std::max(b.t.l,b.ss[j].l);b.y[j].l=std::max(b.y[j].l,b.ss[j].l);b.ss[j].h=std::min({b.ss[j].h,b.y[j].h,b.t.h});if(b.y[j].l>b.y[j].h || b.ss[j].l>b.ss[j].h)return false;}
 for(int j=1;j<H;++j){if(b.y[j-1].l>b.y[j].h)return false;b.y[j-1].h=std::min(b.y[j-1].h,b.y[j].h);b.y[j].l=std::max(b.y[j].l,b.y[j-1].l);}
 return b.t.l<=b.t.h && b.D.l<=b.D.h;
}
I base_at(I t,I mu,I sumy,I sums){
 I sb=mn(t,sums),reff=sumy-mx(ZERO,sums-t),edge=t-sb;
 assert(reff.l>=0 && edge.l>=0);
 I alow=clip(ONE/(mu*KAPPA)-reff,ZERO,edge);
 I ahigh=mx(edge,ONE/(mu*BETA)-reff);
 I low=lg(alow+reff)-mu*KAPPA*alow+mu*BETA*t;
 I high=lg(ahigh+reff)-mu*BETA*ahigh+mu*(t-Z*sb);
 return mx(low,high);
}

// The five functions are individually concave and increasing in (y,w).
// Their maximum is the exact owner net price for this transfer.  We take a
// chord of one function composed with one concave physical score, never of
// the maximum of the scores or of the five functions.
I original_branch_price(int mode,I w,I y,I D){
 I lw=lg(w);
 if(mode==0)return (ONE-ETA)*lw;
 if(mode==1)return lw-lg(ONE+KM*w);
 if(mode==2){assert(y.l>0);return lw+lg(y)-lg(y+KM*D*w);}
 if(mode==3){I x=B*w;return lw+lg(ONE+x)-lg(ONE+x+x*x);}
 assert(mode==4 && y.l>0);
 I v=B*D*w;if(v.h==0)return lw;
 if(y.l>=v.h){I x=v/y;return lw+lg(ONE+x)-lg(ONE+x+x*x);}
 assert(v.l>0);
 I x=y/v;return lw+lg(x)+lg(ONE+x)-lg(ONE+x+x*x);
}
I branch_price(int mode,I w,I y,I D){
 if(mode==5)return rat(499,1000)*localized_g(w);
 return original_branch_price(mode,w,y,D)-HALF_PRICE_CAP;
}
struct Line{I cost,left;};
Line price_line(const Box& box,int j){
 I yl=pt(box.y[j].l),yh=pt(box.y[j].h),t=pt(box.t.l),D=pt(box.D.h),ss=pt(box.ss[j].h);
 int best=0;i64 bestmid=SCALE;
 for(int k=1;k<score_count();++k)if(score_active(k,D,yl)){
  I wl=physical_score(k,t,D,yl,ss),wh=physical_score(k,t,D,yh,ss);
  if(wl.l<=0 || wh.l<=0)continue;
  // Select at the right endpoint so that a score crossing w=1 inside the
  // box is retained even when the midpoint still has zero owner price.
  I wm=wh;if(wm.l>bestmid){bestmid=wm.l;best=k;}
 }
 I wl=pt(physical_score(best,t,D,yl,ss).l),wh=pt(physical_score(best,t,D,yh,ss).l);
 Line ans{ZERO,ZERO};i64 bestvalue=0;
 for(int mode=0;mode<6;++mode){
  if((mode==2 || mode==4) && yl.l==0)continue;
  I vl=branch_price(mode,wl,yl,D),vh=branch_price(mode,wh,yh,D);
  i64 left=std::min(vl.l,vh.l),slope=0;
  if(yh.l>yl.l){i128 raw=i128(vh.l-left)*SCALE/(yh.l-yl.l);slope=i64(std::min(raw,i128(GUARD/4)));}
  I value=pt(left)+pt(slope)*(yh-yl);
  if(value.l>bestvalue){bestvalue=value.l;ans={pt(slope),pt(left)};}
 }
 return ans;
}

long double tangent_guess(const Box& box,const Line lines[4],I mu,bool high,I t){
 long double sums=0,initial=0;for(int j=0;j<H;++j){sums+=fl(pt(box.ss[j].l));initial+=fl(pt(box.y[j].l));}
 long double edge=std::max(0.L,fl(t)-sums),overlap=std::max(0.L,sums-fl(t));initial-=overlap;if(high)initial+=edge;
 struct V{long double cost,width;};std::vector<V> vars;
 vars.push_back({fl(mu)*fl(high?BETA:KAPPA),high?1e100L:edge});
 for(int j=0;j<H;++j){
  long double width=static_cast<long double>(box.y[j].h-box.y[j].l)/SCALE;
  long double cost=fl(lines[j].cost),left=fl(lines[j].left);
  long double zero=(cost>0)?std::max(0.L,std::min(width,-left/cost)):width;
  vars.push_back({0,zero});vars.push_back({cost,width-zero});
 }
 std::sort(vars.begin(),vars.end(),[](V a,V b){return a.cost<b.cost;});
 long double vv=std::max(0.L,initial);
 for(V v:vars){long double w=v.width;if(v.cost>0)w=std::min(w,std::max(0.L,1/v.cost-vv));if(w>0)vv+=w;}
 assert(vv>0 && std::isfinite(vv));return 1/vv;
}
I tangent_bound(const Box& box,const Line lines[4],I mu,bool high,I t,I tangent,bool normalized){
 I sums=ZERO;for(int j=0;j<H;++j)sums=sums+pt(box.ss[j].l);
 I sb=mn(t,sums),overlap=mx(ZERO,sums-t),edge=t-sb;
 I ca=mu*(high?BETA:KAPPA),yy=normalized?tangent*mu:tangent;
 I out=normalized?-lg(tangent):-lg(tangent)-ONE;
 out=out-yy*overlap;
 I aco=yy-ca;
 if(high){if(normalized){assert(tangent.h<=BETA.l);aco.h=std::min<i64>(0,aco.h);}assert(aco.h<=0);out=out+aco*edge;}else out=out+mx(ZERO,aco*edge);
 for(int j=0;j<H;++j){
  I yl=pt(box.y[j].l),yh=pt(box.y[j].h),width=yh-yl;
  I right=lines[j].left+lines[j].cost*width;
  I owner=mx(yy*yl-mx(ZERO,lines[j].left),yy*yh-mx(ZERO,right));
  // P >= max(0,chord).  The additional zero-crossing vertex preserves the
  // genuine w=1 kink instead of forcing tiny boxes around each crossing.
  if(lines[j].cost.l>0 && lines[j].left.l<0 && right.h>=0){
   I crossing=clip(yl-lines[j].left/lines[j].cost,yl,yh);
   owner=mx(owner,yy*crossing);
  }
  out=out+owner;
 }
 if(normalized)out=out-(high?mu*Z*sb:mu*Z*t);
 else out=out+(high?mu*(t-Z*sb):mu*BETA*t);
 return out;
}
i64 rational_tangent(long double value){return std::max<i64>(1,std::min<i64>(GUARD/4,static_cast<i64>(std::min(value*SCALE,static_cast<long double>(GUARD/4)))));}
I ordinary_chord(const Box& box,const Line lines[4],I mu){
 I result=num(-1000),t=pt(box.t.h);
 for(bool high:{false,true}){i64 yy=rational_tangent(tangent_guess(box,lines,mu,high,t));if(high)yy=std::min(yy,(mu*BETA).l);result=mx(result,tangent_bound(box,lines,mu,high,t,pt(yy),false));}
 return result;
}
I normalized_chord(const Box& box,const Line lines[4],I mu){
 // Actual boosted admissibility gives C >= mu*t-1-log(mu)+b.
 // The required +b payment cancels exactly; no uniform b_max loss belongs
 // in this normalized bound. Applying Fenchel to
 // log(mu*(a+reff)) leaves an affine expression in mu.  Evaluate its two
 // endpoints; no unjustified convexity claim about log(mu) is needed.
 I result=num(-1000),mm=pt(mu.l+(mu.h-mu.l)/2),tm=pt(box.t.l+(box.t.h-box.t.l)/2);
 I sums=ZERO;for(int j=0;j<H;++j)sums=sums+pt(box.ss[j].l);
 std::vector<i64> ts{box.t.l,box.t.h,std::min(box.t.h,std::max(box.t.l,sums.l))};
 for(bool high:{false,true}){
  i64 ss=rational_tangent(tangent_guess(box,lines,mm,high,tm)/fl(mm));if(high)ss=std::min(ss,BETA.l);
  I v=num(-1000);
  for(i64 m:{mu.l,mu.h})for(i64 tt:ts)v=mx(v,tangent_bound(box,lines,pt(m),high,pt(tt),pt(ss),true));
  result=mx(result,v);
 }
 return result;
}

struct Plane{I cy,cs,value;};
struct Point{I y,ss;};
I plane_at(const Plane& p,const Box& b,int j,Point x){return p.value+p.cy*(x.y-pt(b.y[j].l))+p.cs*(x.ss-pt(b.ss[j].l));}
I plane_physical_score(int index,I t,I D,I y,I ss){
 if(index<2+int(PROX.size()))return physical_score(index,t,D,y,ss);
 const Planar& p=PLANAR[index-2-PROX.size()];
 // This unenhanced floor is jointly affine in (y,S), hence the score is
 // jointly concave on the whole rectangle without a floor-regime split.
 I total=p.c+p.d,mean=(p.a+p.d)*y+total*t-p.d*ss;
 return mean-gap_upper(total*D,(p.a+p.d)*y);
}
i64 slope_select(i64 value,i64 width){
 if(width==0)return 0;
 i128 q=i128(value)*SCALE/width;
 return i64(std::max(-i128(GUARD/16),std::min(i128(GUARD/16),q)));
}
std::array<Plane,2> price_planes(const Box& b,int j){
 I yl=pt(b.y[j].l),yh=pt(b.y[j].h),sl=pt(b.ss[j].l),sh=pt(b.ss[j].h),t=pt(b.t.l),D=pt(b.D.h);
 std::array<Point,4> vertices{{{yl,sl},{yh,sl},{yl,sh},{yh,sh}}};
 int score=0;i64 best=SCALE;
 for(int k=1;k<score_count();++k)if(score_active(k,D,yl)){
  I wl=plane_physical_score(k,t,D,yl,sh),wr=plane_physical_score(k,t,D,yh,sl);
  if(wl.l>0 && wr.l>best){best=wr.l;score=k;}
 }
 std::array<I,4> ws;
 for(int k=0;k<4;++k)ws[k]=pt(plane_physical_score(score,t,D,vertices[k].y,vertices[k].ss).l);
 std::array<Plane,2> ans;std::array<i64,2> merit{{-GUARD+1,-GUARD+1}};
 for(int mode=0;mode<6;++mode){
  if((mode==2 || mode==4) && yl.l==0)continue;
  std::array<I,4> vals;for(int k=0;k<4;++k)vals[k]=branch_price(mode,ws[k],vertices[k].y,D);
  I cy=pt(slope_select(vals[3].l-vals[2].l,yh.l-yl.l));
  I cs=pt(slope_select(vals[3].l-vals[1].l,sh.l-sl.l));
  I intercept=num(1000);
  for(int k=0;k<4;++k)intercept=mn(intercept,pt(vals[k].l)-cy*(vertices[k].y-yl)-cs*(vertices[k].ss-sl));
  Plane p{cy,cs,pt(intercept.l)};
  int slot=mode==0?0:1;
  if(vals[1].l>merit[slot]){merit[slot]=vals[1].l;ans[slot]=p;}
 }
 return ans;
}
I owner_plane_sup(const Box& b,int j,const Plane& p,I theta,I ell){
 I yl=pt(b.y[j].l),yh=pt(b.y[j].h),sl=pt(b.ss[j].l),sh=pt(b.ss[j].h);
 std::vector<Point> v{{yl,sl},{yh,sl},{yh,sh}};
 if(sh.l>yl.l){v.push_back({sh,sh});v.push_back({yl,yl});}else v.push_back({yl,sh});
 I ans=num(-1000);
 for(Point x:v)ans=mx(ans,theta*x.y-ell*x.ss-mx(ZERO,plane_at(p,b,j,x)));
 for(std::size_t k=0;k<v.size();++k){
  Point a=v[k],z=v[(k+1)%v.size()];I pa=plane_at(p,b,j,a),pz=plane_at(p,b,j,z);
  if(pa.l>0 && pz.l>0)continue;
  if(pa.h<0 && pz.h<0)continue;
  if(pa.l>pz.l){std::swap(a,z);std::swap(pa,pz);}
  I diff=pz-pa;
  if(diff.l<=0){ans=mx(ans,mx(theta*a.y-ell*a.ss,theta*z.y-ell*z.ss));continue;}
  I weight=clip(-pa/diff,ZERO,ONE);
  Point cross{a.y+weight*(z.y-a.y),a.ss+weight*(z.ss-a.ss)};
  ans=mx(ans,theta*cross.y-ell*cross.ss);
 }
 return ans;
}
I normalized_plane_bound(const Box& b,I mu){
 std::array<std::array<Plane,2>,4> planes;
 for(int j=0;j<H;++j)planes[j]=price_planes(b,j);
 // The actual +b payment cancels against boosted admissibility here too.
 // These witnesses are constants, not numerical proof tests.  For each,
 // the true owner supremum is convex in mu; endpoint bounds are valid.
 I best=num(1000);
 for(I s:{pt(BETA.l),pt((BETA*rat(3,4)).l),pt((BETA*rat(1,2)).l)}){
  for(I ell:{mn(pt(Z.l),s),s,pt((s*rat(1,2)).l)}){
   I total=num(-1000);
   for(i64 m:{mu.l,mu.h}){
    I mm=pt(m),theta=s*mm,ll=ell*mm;
    I owner=ZERO;
    for(int j=0;j<H;++j)owner=owner+mn(owner_plane_sup(b,j,planes[j][0],theta,ll),owner_plane_sup(b,j,planes[j][1],theta,ll));
    I g=mx(ell-Z,s-BETA);
    I value=-lg(s)+mx(mm*pt(b.t.l)*g,mm*pt(b.t.h)*g)+owner;
    total=mx(total,value);
   }
   best=mn(best,total);
  }
 }
 return best;
}

std::uint64_t BOUND_CALLS=0;
i64 upper_bound(const Box& box,TDInfo td){
 ++BOUND_CALLS;
 I sumy=ZERO,sums=ZERO,prices=ZERO,q=ZERO;
 for(int j=0;j<H;++j){sumy=sumy+pt(box.y[j].h);sums=sums+pt(box.ss[j].l);q=mx(q,pt(box.ss[j].l)-BETA*pt(box.y[j].h));I w=physical(pt(box.t.l),pt(box.D.h),pt(box.y[j].l),pt(box.ss[j].h));prices=prices+owner_net(pt(w.l),pt(box.y[j].l),pt(box.D.h));}
 q=mn(q,Z*pt(box.t.l));q=pt(std::max<i64>(0,q.l));
 I direct=num(-1000);
 std::array<I,2> credits;
 int n=0;for(i64 m:{td.mu.l,td.mu.h}){I mu=pt(m),credit=home_C(pt(box.t.l),pt(box.D.l),q,mu);credits[n++]=credit;direct=mx(direct,base_at(pt(box.t.h),mu,sumy,sums)-prices-credit);}
 direct=direct-LAM*num(H-1)+PRIVATE_PAYMENT_CAP;if(direct.h<U.l)return direct.h;
 Line lines[4];for(int j=0;j<H;++j)lines[j]=price_line(box,j);
 I chord=num(-1000);n=0;for(i64 m:{td.mu.l,td.mu.h})chord=mx(chord,ordinary_chord(box,lines,pt(m))-credits[n++]);
 chord=chord-LAM*num(H-1)+PRIVATE_PAYMENT_CAP;
 i64 ans=std::min(direct.h,chord.h);
 if(ans>=U.l){I normalized=normalized_chord(box,lines,td.mu)-LAM*num(H-1);ans=std::min(ans,normalized.h);}
  if(ans>=U.l){I planar=normalized_plane_bound(box,td.mu)-LAM*num(H-1);ans=std::min(ans,planar.h);}
 return ans;
}

I& coordinate(Box& b,int d){if(d==0)return b.t;if(d==1)return b.D;return d%2?b.ss[(d-2)/2]:b.y[(d-2)/2];}
int choose_dimension(const Box& b){
 long double best=-1;int dim=-1;
 auto take=[&](I v,long double factor,int k){long double width=static_cast<long double>(v.h-v.l)/SCALE;long double score=width*factor;if(score>best){best=score;dim=k;}};
 take(b.t,4,0);take(b.D,1/(1+fl(b.D)),1);
 for(int j=0;j<H;++j){take(b.y[j],1.7L/(1+fl(b.y[j])),2+2*j);take(b.ss[j],1.1L,3+2*j);}
 assert(dim>=0 && best>0);return dim;
}
void print_box(const Box& b){std::cout<<"BOX depth "<<b.depth<<" t "<<b.t.l<<" "<<b.t.h<<" D "<<b.D.l<<" "<<b.D.h;for(int j=0;j<H;++j)std::cout<<" y"<<j<<" "<<b.y[j].l<<" "<<b.y[j].h<<" S"<<j<<" "<<b.ss[j].l<<" "<<b.ss[j].h;std::cout<<"\n";}
int main(int argc,char** argv){
 if(argc<5){std::cerr<<"usage: verify_B_simple tables.txt input.txt h maxboxes [tlo thi Dlo Dhi ymax] (bounds in decimal, converted outward)\n";return 2;}
 read_tables(argv[1]);read_input(argv[2]);assert(PROX.size()==5 && PLANAR.size()==1 && score_count()==8);
 if(std::string(argv[3])=="evaluate"){
  assert(argc==5);std::ifstream tests(argv[4]);assert(tests);std::string kind;i64 a,b,c,d;
  while(tests>>kind>>a>>b>>c>>d){I value=kind=="C"?home_C(pt(a),pt(b),pt(c),pt(d)):owner_net(pt(a),pt(b),pt(c));std::cout<<kind<<" "<<value.l<<" "<<value.h<<"\n";}
  return 0;
 }
 if(std::string(argv[3])=="multiplier"){
  assert(argc==8);I tt(std::stoll(argv[4]),std::stoll(argv[5])),dd(std::stoll(argv[6]),std::stoll(argv[7]));
  TDInfo v=td_info(tt,dd);std::cout<<"MULTIPLIER "<<v.mu.l<<" "<<v.mu.h<<" "<<v.excluded<<"\n";return 0;
 }
 if(std::string(argv[3])=="credit"){
  assert(argc==8);I value=home_C(pt(std::stoll(argv[4])),pt(std::stoll(argv[5])),pt(std::stoll(argv[6])),pt(std::stoll(argv[7])));
  std::cout<<"CREDIT "<<value.l<<" "<<value.h<<"\n";return 0;
 }
 if(std::string(argv[3])=="price"){
  assert(argc==7);I value=owner_net(pt(std::stoll(argv[4])),pt(std::stoll(argv[5])),pt(std::stoll(argv[6])));
  std::cout<<"PRICE "<<value.l<<" "<<value.h<<"\n";return 0;
 }
 H=std::stoi(argv[3]);assert(1<=H && H<=4);std::uint64_t maxboxes=std::stoull(argv[4]);
 Box initial;initial.t=I(0,TMAX.h);initial.D=I(0,DMAX.h);I ymax=num(128);
 auto decimal_bound=[](const char* s){std::string x(s);auto pos=x.find('.');i64 n=0,d=1;for(char c:x)if(c!='.'){assert(c>='0' && c<='9');n=n*10+(c-'0');}if(pos!=std::string::npos)for(std::size_t j=pos+1;j<x.size();++j)d*=10;return rat(n,d);};
 if(argc>=9){initial.t=I(decimal_bound(argv[5]).l,decimal_bound(argv[6]).h);initial.D=I(decimal_bound(argv[7]).l,decimal_bound(argv[8]).h);}if(argc>=10)ymax=decimal_bound(argv[9]);
 for(int j=0;j<H;++j){initial.y[j]=I(0,ymax.h);initial.ss[j]=I(0,initial.t.h);}
 // The exact A rectangle is t<=3/25 and D<=(21/20)*R4^2.
 // Two outward-rounded roots cover its complement. Their tiny boundary
 // overlap is intentional; no unproved point is excluded by rounding.
 I A_T=rat(3,25),R4=rat(100032357,100000000),A_D=rat(21,20)*R4*R4;
 std::vector<Box> stack;
 if(initial.D.h>=A_D.l){Box highD=initial;highD.D.l=std::max(highD.D.l,A_D.l);highD.t.h=std::min(highD.t.h,A_T.h);if(highD.t.l<=highD.t.h)stack.push_back(highD);}
 if(initial.t.h>=A_T.l){Box highT=initial;highT.t.l=std::max(highT.t.l,A_T.l);if(highT.t.l<=highT.t.h)stack.push_back(highT);}
 const std::uint64_t roots=stack.size();
 std::uint64_t count=0,accepted=0,excluded=0,analytic_band=0;int maxdepth=0;i64 highest=-GUARD+1;
 auto start=std::chrono::steady_clock::now();
 std::cout<<"START B_SIMPLE h "<<H<<" scale "<<SCALE<<" U "<<U.l<<" "<<U.h<<" tmax "<<TMAX.h<<" Dmax "<<DMAX.h<<" roots "<<roots<<" A_D "<<A_D.l<<" "<<A_D.h<<"\n";print_box(initial);std::cout.flush();
 while(!stack.empty()){
  Box b=stack.back();stack.pop_back();++count;maxdepth=std::max(maxdepth,b.depth);
  if(count%500000==0){double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"PROGRESS "<<count<<" pending "<<stack.size()<<" depth "<<maxdepth<<" analytic_band "<<analytic_band<<" seconds "<<elapsed<<"\n";std::cout.flush();}
  if(!propagate(b)){++excluded;continue;}
  TDInfo td=td_info(b.t,b.D);if(td.excluded){++excluded;continue;}
  // One standalone theorem for every H=1,2,3,4; all tests point inward.
  // It proves the SAME boosted target and capped-price debit with reserve 1e-6.
  if(b.t.l>=rat(1199,10000).h && b.t.h<=rat(3001,10000).l && b.D.h<=rat(23,20).l){++accepted;++analytic_band;continue;}
  // Every other accepted B box has its own directed scalar bound.
  i64 bound=upper_bound(b,td);
  if(bound<U.l){++accepted;highest=std::max(highest,bound);continue;}
  if(count>=maxboxes || b.depth>=180){std::cout<<"INCOMPLETE boxes "<<count<<" pending "<<stack.size()+1<<" accepted "<<accepted<<" excluded "<<excluded<<" bound "<<bound<<" mu "<<td.mu.l<<" "<<td.mu.h<<"\n";print_box(b);return 1;}
  int d=choose_dimension(b);I v=coordinate(b,d);i64 mid=v.l+(v.h-v.l)/2;assert(mid>v.l && mid<v.h);Box left=b,right=b;coordinate(left,d).h=mid;coordinate(right,d).l=mid;left.depth=right.depth=b.depth+1;stack.push_back(right);stack.push_back(left);
 }
 double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
 std::cout<<"COMPLETE B_SIMPLE h "<<H<<" boxes "<<count<<" accepted "<<accepted<<" excluded "<<excluded<<" analytic_band "<<analytic_band<<" roots "<<roots<<" depth "<<maxdepth<<" highest_upper "<<highest<<" / "<<SCALE<<" threshold_lower "<<U.l<<" / "<<SCALE<<" seconds "<<elapsed<<"\n";
 return 0;
}
