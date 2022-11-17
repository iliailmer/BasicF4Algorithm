BasicF4Alg := proc(Ff, vars)
    local Gg, ordering_, s, Bb, Bb_prime, i, j,
    lcms, idxs, LM, LM_i, LM_j, Mm, LCM;

    Gg := Ff;
    ordering_ := tdeg(op(vars)):
    s := numelems(Ff); # s in the Cox, Little, O'Shea book
    t := s;
    Bb := [];
    for j from 1 to s do
        for i from 1 to j do
            if i<j then
                Bb := [op(Bb), [i,j]];
            end if;
        end do;
    end do;

    while numelems(Bb)>0 do
        lcms := []:
        # pick pairs that have minimal degree in Ff
        Bb_prime := Selection(Bb, Gg, ordering_);
        Bb := select(x->not (x in Bb_prime), Bb);
        for idxs in Bb_prime do
            i := idxs[1];
            j := idxs[2];
            LM_i := Groebner:-LeadingMonomial(Gg[i], ordering_);
            LM_j := Groebner:-LeadingMonomial(Gg[j], ordering_);
            LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Gg[i], ordering_)[2]*Gg[i];
            lcms := [op(lcms), LCM];
            LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Gg[j], ordering_)[2]*Gg[j];
            lcms := [op(lcms), LCM];
        end do;
        lcms := [op({op(lcms)})];
        Mm, Mon_, Hh := ComputeM(lcms, Gg, ordering_):
        N := LinearAlgebra:-ReducedRowEchelonForm(Mm);
        leading_monomials_ := {seq(Groebner:-LeadingMonomial(each, ordering_),  each in Hh)};

        N := convert(N . Vector(Mon_), list);
        N_plus := [];
        for each in N do
            divisibility:= [seq(divide(Groebner:-LeadingMonomial(each, ordering_), lm), lm in leading_monomials_)]:
            if not (true in divisibility) then
                N_plus := [op(N_plus), each]:
            end if:
        end do:
        N_plus := select(x->not type(x, numeric), N_plus);
        for each in N_plus do
            Gg := [op(Gg), each]:
            t:= t+1;
            for i from 1 to t-1 do
                Bb := [op(Bb), [i,t]]:
            end do:
        end do:
    end do:
    return Gg;
end proc:

GetMonomials := proc(Hh, ordering_)
    local Mon_, f, each:
    Mon_ := [seq(op(expand(f)), f in Hh)];
    Mon_ := [seq(each/Groebner:-LeadingCoefficient(each, ordering_), each in Mon_)];
    Mon_ := [op(sort(add(Mon_), order=ordering_))]:
    Mon_ := [seq(each/Groebner:-LeadingCoefficient(each, ordering_), each in Mon_)];
    return Mon_
end proc:

ComputeM := proc(lcms, Gg, ordering_)
    local Hh, done_, Mon_, tmp, x_beta, idx, each, Mm, broken_poly, r, c, f;
    Hh := lcms;
    done_ := {seq(Groebner:-LeadingMonomial(each, ordering_), each in Hh)};
    Mon_ := GetMonomials(Hh, ordering_);
    while ({op(done_)} <> {op(Mon_)}) do
        tmp := select(x->not x in done_, Mon_);
        if numelems(tmp) = 0 then
            break:
        end if:
        x_beta := tmp[1];
        done_ := [op(done_), x_beta];
        for each in Gg do
            if divide(x_beta, Groebner:-LeadingMonomial(each, ordering_)) then
                Hh := [op(Hh), each * x_beta/Groebner:-LeadingMonomial(each, ordering_)];
                break;
            end if;
        end do;
        Mon_ := GetMonomials(Hh, ordering_);
    od:
    # construct a matrix of coefficients of monomials in Hh
    Mm := Matrix(numelems(Hh), numelems(Mon_)):

    for r from 1 to numelems(Hh) do
        for c from 1 to numelems(Mon_) do
            broken_poly := [op(expand(Hh[r]))];
            for each in broken_poly do
                if type(each/Mon_[c], numeric) then
                    Mm[r,c] := each/Mon_[c];
                    break;
                else
                    Mm[r,c] := 0;
                end if;
            end do;
        end do;
    end do;

    return Mm, Mon_, Hh:
end proc:

Selection := proc(Bb, Ff, ordering_)
    local degrees, dd, p;
    degrees := [
        seq(
            degree(
                lcm(Groebner:-LeadingMonomial(Ff[p[1]], ordering_), Groebner:-LeadingMonomial(Ff[p[2]], ordering_))
                ),
        p in Bb)
    ]:
    dd := min(degrees):
    return select(x->degree(lcm(Groebner:-LeadingMonomial(Ff[x[1]], ordering_), Groebner:-LeadingMonomial(Ff[x[2]], ordering_)))=dd, Bb):
end proc:
