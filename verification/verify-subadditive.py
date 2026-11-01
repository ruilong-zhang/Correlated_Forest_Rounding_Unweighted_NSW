#!/usr/bin/env  python3
"""Exact rational verification for the subadditive concentration and calibration appendix.
Run  with  Python  3.10  or  later.  No  third-party  packages  are  needed.
Every  decisive  comparison  uses  Fraction  arithmetic.
"""
from  fractions  import  Fraction  as  Q
from  functools  import  lru_cache
from  math  import  factorial,  isqrt
import  json
from  pathlib  import  Path
TERMS  =  18
SCALE  =  1  <<  54
def  ceilq(x):
        return  -((-x.numerator)  //  x.denominator)
def  outward(lo,  hi):
        return  (Q((lo  *  SCALE).numerator  //  (lo  *  SCALE).denominator,
                            SCALE),  Q(ceilq(hi  *  SCALE),  SCALE))
def  unit_log(x):
        assert  1  <=  x  <=  2
        z  =  (x  -  1)  /  (x  +  1)
        z2,  power,  total  =  z  *  z,  z,  Q(0)
        for  j  in  range(TERMS):
                total  +=  2  *  power  /  (2  *  j  +  1)
                power  *=  z2
        tail  =  2  *  power  /  ((2  *  TERMS  +  1)  *  (1  -  z2))
        return  total,  total  +  tail
L2,  U2  =  unit_log(Q(2))
@lru_cache(None)
def  log_interval(x):
        assert  x  >  0
        k  =  x.numerator.bit_length()  -  x.denominator.bit_length()
        power  =  Q(2  **  k)  if  k  >=  0  else  Q(1,  2  **  (-k))
        t  =  x  /  power
        if  t  <  1:
                k,  t  =  k  -  1,  2  *  t
        elif  t  >=  2:
                k,  t  =  k  +  1,  t  /  2
        lo,  hi  =  unit_log(t)
        if  k  >=  0:
                lo,  hi  =  lo  +  k  *  L2,  hi  +  k  *  U2
        else:
                lo,  hi  =  lo  +  k  *  U2,  hi  +  k  *  L2
        return  outward(lo,  hi)
def  exp_interval(x,  terms=48):
        assert  0  <=  x  <  terms  +  2
        term  =  total  =  Q(1)
        for  k  in  range(1,  terms  +  1):
                term  *=  x  /  k
                total  +=  term
        first  =  term  *  x  /  (terms  +  1)
        return  outward(total,  total  +  first  /  (1  -  x  /  (terms  +  2)))
def  e_interval():
        total  =  sum((Q(1,  factorial(k))  for  k  in  range(37)),  Q(0))
        tail  =  Q(1,  factorial(37))  /  (1  -  Q(1,  38))
        return  total,  total  +  tail
def  sqrt_interval(x):
        d  =  1  <<  72
        a  =  x  *  d  *  d
        k  =  isqrt(a.numerator  //  a.denominator)
        lo,  hi  =  Q(k,  d),  Q(k  +  1,  d)
        assert  lo  *  lo  <=  x  <=  hi  *  hi
        return  lo,  hi
def  concentration():
        #  Upper  enclosures  for  b(j/100)  in  the  quantile  bound.
        grid  =  []
        for  j  in  range(1,  100):
                u  =  Q(j,  100)
                upper  =  log_interval((1  +  u  +  u  *  u)  /  (u  *  u))[1]
                b  =  (1  -  u)  *  upper  /  L2
                grid.append((u,  Q(ceilq(b  *  SCALE),  SCALE)))
        records  =  []
        z  =  2375
        while  z  <  40000:
                nxt  =  (101  *  z  +  99)  //  100
                lower,  upper  =  Q(z,  1000),  Q(nxt,  1000)
                value  =  Q(0)
                for  u,  b  in  grid:
                        t  =  (lower  -  b)  /  (2  -  u)
                        if  t  >  1:
                                value  +=  log_interval(t)[0]  /  100
                needed  =  log_interval(upper  /  Q(19,  8))[1]
                margin  =  value  -  needed
                assert  margin  >  0
                records.append({"left":  str(lower),  "right":  str(upper),
                                                "margin":  str(margin)})
                z  =  nxt
        tail  =  Q(3,  1024)  +  Q(3,  20)
        assert  tail  <  Q(87,  512)  <  log_interval(Q(19,  16))[0]
        assert  L2  >  Q(2,  3)  and  log_interval(Q(10))[1]  <  3
        assert  len(records)  ==  282
        smallest  =  min(records,  key=lambda  a:  Q(a["margin"]))
        return  records,  smallest
def  approximation_constants():
        elo,  ehi  =  e_interval()
        l2,  u2  =  log_interval(Q(2))
        clo,  chi  =  2  /  (ehi  *  u2),  2  /  (elo  *  l2)
        lower  =  upper  =  Q(0)
        for  k  in  range(31):
                lo,  hi  =  log_interval(Q(k  +  1))
                lower  +=  lo  /  factorial(k)
                upper  +=  hi  /  factorial(k)
        lamlo  =  lower  /  ehi
        lamhi  =  upper  /  elo  +  Q(31,  30)  /  (elo  *  factorial(30))
        alo  =  clo  *  exp_interval(lamlo)[0]
        ahi  =  chi  *  exp_interval(lamhi)[1]
        rho,  s  =  Q(4,  15),  Q(11,  15)
        corelo  =  alo  /  rho  *  exp_interval(s  /  (2  *  ahi))[0]
        corehi  =  ahi  /  rho  *  exp_interval(s  /  (2  *  alo))[1]
        rtlo,  rthi  =  sqrt_interval(Q(30))
        total_lo  =  1  +  corelo  +  Q(19,  4)  *  rtlo
        total_hi  =  1  +  corehi  +  Q(19,  4)  *  rthi
        assert  total_hi  <  Q(178,  5)
        #  Rational  bounds  quoted  in  the  text.
        assert  Q(106147,  100000)  <  clo  <  chi  <  Q(106148,  100000)
        assert  Q(57340,  100000)  <  lamlo  <  lamhi  <  Q(57341,  100000)
        assert  Q(3559743,  100000)  <  total_lo
        assert  total_hi  <  Q(3559744,  100000)
        return  {"C":  [str(clo),  str(chi)],
                        "Lambda":  [str(lamlo),  str(lamhi)],
                        "subadditive":  [str(total_lo),  str(total_hi)],
                        "margin_below_35_6":  str(Q(178,  5)  -  total_hi)}
def  main():
        records,  smallest  =  concentration()
        constants  =  approximation_constants()
        result  =  {"status":  "PASS",  "interval_count":  len(records),
                            "minimum":  smallest,  "constants":  constants,
                            "intervals":  records}
        Path("certificate.json").write_text(json.dumps(result,  indent=2)  +  "\n")
        print("PASS:  282  concentration  intervals  and  all  constant  bounds")
        print("Minimum  concentration  margin:",  smallest["margin"])
        print("Interval:",  smallest["left"],  smallest["right"])
        print("Exact  output  written  to  certificate.json")
if  __name__  ==  "__main__":
        main()
