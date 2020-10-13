#	here are the guts of the fast extended products

sta[geometric_product] := proc(a::multivector,b::multivector)::multivector;
        local a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16\
             ,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,sigma;
        dirac[sparts](a, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16);
        dirac[sparts](b, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

        sigma:=[a1*b1+a2*b2-a3*b3-a4*b4-a5*b5-a6*b6-a7*b7-a8*b8+a9*b9+a10*b10+a11*b11+a12*b12-a13*b13-a14*b14-a15*b15-a16*b16, a1*b2
+a2*b1-a3*b9-a4*b10-a5*b11-a6*b13-a7*b14-a8*b15+a9*b3+a10*b4+a11*b5+a12*b16-a13*b6-a14*b7-a15*b8-a16*b12, a1*b3-a2*b9
+a3*b1-a4*b8+a5*b7-a6*b12-a7*b5+a8*b4+a9*b2-a10*b15+a11*b14-a12*b6+a13*b16+a14*b11-a15*b10-a16*b13, a1*b4-a2*b10+a3*
b8+a4*b1-a5*b6+a6*b5-a7*b12-a8*b3+a9*b15+a10*b2-a11*b13-a12*b7-a13*b11+a14*b16+a15*b9-a16*b14, a1*b5-a2*b11-a3*b7+a4*
b6+a5*b1-a6*b4+a7*b3-a8*b12-a9*b14+a10*b13+a11*b2-a12*b8+a13*b10-a14*b9+a15*b16-a16*b15, a1*b6+a2*b13-a3*b12-a4*b5+a5
*b4+a6*b1-a7*b8+a8*b7+a9*b16+a10*b11-a11*b10-a12*b3+a13*b2-a14*b15+a15*b14+a16*b9, a1*b7+a2*b14+a3*b5-a4*b12-a5*b3+a6
*b8+a7*b1-a8*b6-a9*b11+a10*b16+a11*b9-a12*b4+a13*b15+a14*b2-a15*b13+a16*b10, a1*b8+a2*b15-a3*b4+a4*b3-a5*b12-a6*b7+a7
*b6+a8*b1+a9*b10-a10*b9+a11*b16-a12*b5-a13*b14+a14*b13+a15*b2+a16*b11, a1*b9-a2*b3+a3*b2-a4*b15+a5*b14-a6*b16-a7*b11+
a8*b10+a9*b1-a10*b8+a11*b7-a12*b13+a13*b12+a14*b5-a15*b4-a16*b6, a1*b10-a2*b4+a3*b15+a4*b2-a5*b13+a6*b11-a7*b16-a8*b9
+a9*b8+a10*b1-a11*b6-a12*b14-a13*b5+a14*b12+a15*b3-a16*b7, a1*b11-a2*b5-a3*b14+a4*b13+a5*b2-a6*b10+a7*b9-a8*b16-a9*b7
+a10*b6+a11*b1-a12*b15+a13*b4-a14*b3+a15*b12-a16*b8, a1*b12-a2*b16+a3*b6+a4*b7+a5*b8+a6*b3+a7*b4+a8*b5+a9*b13+a10*b14
+a11*b15+a12*b1-a13*b9-a14*b10-a15*b11+a16*b2, a1*b13+a2*b6-a3*b16-a4*b11+a5*b10+a6*b2-a7*b15+a8*b14+a9*b12+a10*b5-
a11*b4-a12*b9+a13*b1-a14*b8+a15*b7+a16*b3, a1*b14+a2*b7+a3*b11-a4*b16-a5*b9+a6*b15+a7*b2-a8*b13-a9*b5+a10*b12+a11*b3-
a12*b10+a13*b8+a14*b1-a15*b6+a16*b4, a1*b15+a2*b8-a3*b10+a4*b9-a5*b16-a6*b14+a7*b13+a8*b2+a9*b4-a10*b3+a11*b12-a12*
b11-a13*b7+a14*b6+a15*b1+a16*b5, a1*b16-a2*b12+a3*b13+a4*b14+a5*b15+a6*b9+a7*b10+a8*b11+a9*b6+a10*b7+a11*b8+a12*b2-
a13*b3-a14*b4-a15*b5+a16*b1]
end proc;

sta[interior_product] := proc(a::multivector,b::multivector)::multivector;
        local a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16\
             ,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,sigma;
        dirac[sparts](a, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16);
        dirac[sparts](b, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

        sigma:=[a1*b1+a2*b2-a3*b3-a4*b4-a5*b5-a6*b6-a7*b7-a8*b8+a9*b9+a10*b10+a11*b11+a12*b12-a13*b13-a14*b14-a15*b15-a16*b16, a1*b2
+a2*b1-a3*b9-a4*b10-a5*b11-a6*b13-a7*b14-a8*b15+a9*b3+a10*b4+a11*b5+a12*b16-a13*b6-a14*b7-a15*b8-a16*b12, a1*b3-a2*b9
+a3*b1-a4*b8+a5*b7-a6*b12-a7*b5+a8*b4+a9*b2-a10*b15+a11*b14-a12*b6+a13*b16+a14*b11-a15*b10-a16*b13, a1*b4-a2*b10+a3*
b8+a4*b1-a5*b6+a6*b5-a7*b12-a8*b3+a9*b15+a10*b2-a11*b13-a12*b7-a13*b11+a14*b16+a15*b9-a16*b14, a1*b5-a2*b11-a3*b7+a4*
b6+a5*b1-a6*b4+a7*b3-a8*b12-a9*b14+a10*b13+a11*b2-a12*b8+a13*b10-a14*b9+a15*b16-a16*b15, a1*b6-a12*b3+a13*b2+a16*b9+
a2*b13-a3*b12+a6*b1+a9*b16, a1*b7+a10*b16-a12*b4+a14*b2+a16*b10+a2*b14-a4*b12+a7*b1, a1*b8+a11*b16-a12*b5+a15*b2+a16*
b11+a2*b15-a5*b12+a8*b1, a1*b9+a14*b5-a15*b4-a16*b6-a4*b15+a5*b14-a6*b16+a9*b1, a1*b10+a10*b1-a13*b5+a15*b3-a16*b7+a3
*b15-a5*b13-a7*b16, a1*b11+a11*b1+a13*b4-a14*b3-a16*b8-a3*b14+a4*b13-a8*b16, a1*b12+a12*b1+a16*b2-a2*b16, a1*b13+a13*
b1+a16*b3-a3*b16, a1*b14+a14*b1+a16*b4-a4*b16, a1*b15+a15*b1+a16*b5-a5*b16, a1*b16+a16*b1]
end proc;

sta[exterior_product] := proc(a::multivector,b::multivector)::multivector;
        local a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16\
             ,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,sigma;
        dirac[sparts](a, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16);
        dirac[sparts](b, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

        sigma:=[a1*b1, a1*b2+a2*b1, a1*b3+a3*b1, a1*b4+a4*b1, a1*b5+a5*b1, a1*b6-a4*b5+a5*b4+a6*b1, a1*b7+a3*b5-a5*b3+a7*b1, a1*b8-
a3*b4+a4*b3+a8*b1, a1*b9-a2*b3+a3*b2+a9*b1, a1*b10+a10*b1-a2*b4+a4*b2, a1*b11+a11*b1-a2*b5+a5*b2, a1*b12+a12*b1+a3*b6
+a4*b7+a5*b8+a6*b3+a7*b4+a8*b5, a1*b13+a10*b5-a11*b4+a13*b1+a2*b6-a4*b11+a5*b10+a6*b2, a1*b14+a11*b3+a14*b1+a2*b7+a3*
b11-a5*b9+a7*b2-a9*b5, a1*b15-a10*b3+a15*b1+a2*b8-a3*b10+a4*b9+a8*b2+a9*b4, a16*b1+a1*b16+a10*b7+a11*b8+a12*b2-a13*b3
-a14*b4-a15*b5-a2*b12+a3*b13+a4*b14+a5*b15+a6*b9+a7*b10+a8*b11+a9*b6]
end proc;

sta[commutator_product] := proc(a::multivector,b::multivector)::multivector;
        local a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16\
             ,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,sigma;
        dirac[sparts](a, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16);
        dirac[sparts](b, b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16);

        sigma:=[0, a10*b4+a11*b5+a12*b16-a16*b12-a3*b9-a4*b10-a5*b11+a9*b3, a13*b16-a16*b13-a2*b9-a4*b8+a5*b7-a7*b5+a8*b4+a9*b2, a10
*b2+a14*b16-a16*b14-a2*b10+a3*b8-a5*b6+a6*b5-a8*b3, a11*b2+a15*b16-a16*b15-a2*b11-a3*b7+a4*b6-a6*b4+a7*b3, a10*b11-
a11*b10-a14*b15+a15*b14-a4*b5+a5*b4-a7*b8+a8*b7, a11*b9+a13*b15-a15*b13+a3*b5-a5*b3+a6*b8-a8*b6-a9*b11, -a10*b9-a13*
b14+a14*b13-a3*b4+a4*b3-a6*b7+a7*b6+a9*b10, -a10*b8+a11*b7-a12*b13+a13*b12-a2*b3+a3*b2-a7*b11+a8*b10, -a11*b6-a12*b14
+a14*b12-a2*b4+a4*b2+a6*b11-a8*b9+a9*b8, a10*b6-a12*b15+a15*b12-a2*b5+a5*b2-a6*b10+a7*b9-a9*b7, a10*b14+a11*b15-a13*
b9-a14*b10-a15*b11+a16*b2-a2*b16+a9*b13, -a12*b9-a14*b8+a15*b7+a16*b3-a3*b16-a7*b15+a8*b14+a9*b12, a10*b12-a12*b10+
a13*b8-a15*b6+a16*b4-a4*b16+a6*b15-a8*b13, a11*b12-a12*b11-a13*b7+a14*b6+a16*b5-a5*b16-a6*b14+a7*b13, a12*b2-a13*b3-
a14*b4-a15*b5-a2*b12+a3*b13+a4*b14+a5*b15]
end proc;
