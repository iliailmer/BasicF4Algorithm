BasicF4Alg := proc(Ff, vars, {weights:="", max_iter:=-1})
    local Gg, ordering_, s, Bb, Bb_prime, i, j, 
    lcms, idxs, LM, LM_i, LM_j, Mm, LCM, numiter,
    global_counts, v, t, start, dd, pair, each, Mon_, Hh, N,
    N_plus, leading_monomials_, divisibility, lm, all_vars:

    global_counts := table([seq(v=0, v in vars)]):
    Gg := Ff:
    ordering_ := tdeg(op(vars)):
    s := numelems(Ff): # s in the Cox, Little, O'Shea book
    t := s:
    Bb := []:
    for j from 1 to s do
        for i from 1 to j do
            if i<j then
                Bb := [op(Bb), [i,j]]:
            end if:
        end do:
    end do:
    numiter := 0:
    while numelems(Bb)>0 do
        printf(`\n\niteration no. %a\t max. iterations: %a\n`, numiter+1, max_iter):
        numiter := numiter + 1:
        lcms := []:
        # pick pairs that have minimal degree in Ff
        start := time():
        dd, Bb_prime := Selection(Bb, Gg, ordering_):
        Bb := select(x->not (x in Bb_prime), Bb):
        for idxs in Bb_prime do
            i := idxs[1]:
            j := idxs[2]:
            LM_i := Groebner:-LeadingMonomial(Gg[i], ordering_):
            LM_j := Groebner:-LeadingMonomial(Gg[j], ordering_):
            LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Gg[i], ordering_)[2]*Gg[i]:
            lcms := [op(lcms), LCM]:
            LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Gg[j], ordering_)[2]*Gg[j]:
            lcms := [op(lcms), LCM]:
        end do:
        printf("time for selection: %a\n", time()-start):
        lcms := [op({op(lcms)})]:
        for pair in GetVarCount({seq(Groebner:-LeadingMonomial(each, ordering_), each in lcms)}) do
            global_counts[lhs(pair)] := global_counts[lhs(pair)] + rhs(pair):
        end do;
        if numelems(select(x->rhs(x)>0, [entries(global_counts,`pairs`)]))>1*numelems(global_counts)/2 then
            return global_counts:
        end if:
        start := time():
        Mm, Mon_, Hh := ComputeM(lcms, Gg, ordering_):
        writeto(sprintf(`./Matrices/M%s_iteration_%a.mpl`, weights, numiter)):
            printf("read \"/home/iilmer/SIAN-wdeg-new/structure_search/systems/regular_sequence_check/custom_f4.mpl\":\n"):
            printf("read \"/home/iilmer/SIAN-wdeg-new/structure_search/systems/regular_sequence_check/reg_seq.mpl\":\n"):
            printf(`dd:= %a:\n`, dd):
            printf(`done_before := %a:\n`, {seq(Groebner:-LeadingMonomial(each, ordering_), each in lcms)}):
            printf(`regseq := RegSeq([op(done_before)], [], 1):\n`):
            printf(`Mm := %a:\n`, Mm):
            printf(`#done_after\nMon_ := %a:\n`, Mon_):
            printf(`Hh := %a:\n`, Hh):
            # printf(`regseq := numelems(RegSeq(Mon_, [], 1)):\n`):
        writeto(terminal):
        printf("time for computing M: %a, size:%a\n", time()-start, ArrayTools[Size](Mm)):
        start := time():
        N := LinearAlgebra:-ReducedRowEchelonForm(Mm):
        printf("time for computing N (rref): %a\n", time()-start):
        start := time():
        leading_monomials_ := {seq(Groebner:-LeadingMonomial(each, ordering_),  each in Hh)}:
        printf("time for computing leading monomials: %a\n", time()-start):

        N := convert(N . Vector(Mon_), list):
        N_plus := []:
        start := time():
        for each in N do
            divisibility:= [seq(divide(Groebner:-LeadingMonomial(each, ordering_), lm), lm in leading_monomials_)]:
            if not (true in divisibility) then
                N_plus := [op(N_plus), each]:
            end if:
        end do:
        printf("time for computing N+: %a\n", time()-start):
        N_plus := select(x->not type(x, numeric), N_plus):
        for each in N_plus do
            Gg := [op(Gg), each]:
            t:= t+1:
            for i from 1 to t-1 do
                Bb := [op(Bb), [i,t]]:
            end do:
        end do:
        appendto(sprintf(`./Matrices/M%s_iteration_%a.mpl`, weights, numiter)):
            printf(`Gg := %a:\n`, Gg):
        writeto(terminal):
        if numiter = max_iter then
            # print("HERE"):
            break:
        end if:
    end do:
    return Gg:
end proc:

GetVarCount := proc(monomials, global_counts_table)
    local all_vars, counts, var, mon, v;
    all_vars := indets(monomials):
    counts := table([seq(v=0, v in all_vars)]):
    for v in all_vars do
        for mon in monomials do
            if v in {op(mon)} then
            counts[v] := counts[v]+1:
            end if:
        end do:
    end do:
    return [entries(counts, `pairs`)];
end proc:

GetMonomials := proc(Hh, ordering_)
    local Mon_, f, each:
    Mon_ := [seq(op(expand(f)), f in Hh)]:
    Mon_ := [seq(each/Groebner:-LeadingCoefficient(each, ordering_), each in Mon_)]:
    Mon_ := [op(sort(add(Mon_), order=ordering_))]:
    Mon_ := [seq(each/Groebner:-LeadingCoefficient(each, ordering_), each in Mon_)]:
    return Mon_
end proc:

ComputeM := proc(lcms, Gg, ordering_)
    local Hh, done_, Mon_, tmp, x_beta, idx, each, Mm, broken_poly, r, c, f, numiter, start:
    numiter := 0:
    Hh := lcms: # comonents of S-polynomials
    done_ := {seq(Groebner:-LeadingMonomial(each, ordering_), each in Hh)}: # collection of leading monomials
    Mon_ := GetMonomials(Hh, ordering_): # all monomials
    
    while ({op(done_)} <> {op(Mon_)}) do
        # printf(`\t\tComputeM iteration no. %a\n`, numiter+1):
        numiter := numiter + 1:
        tmp := select(x->not x in done_, Mon_): # select monomials that are not in leading monomial collection
        # print({op(Mon_)} minus {op(done_)}):
        if numelems(tmp) = 0 then
            break:
        end if:
        x_beta := tmp[1]:
        done_ := [op(done_), x_beta]: # add first such monomial ordered by our ordering
        start := time():
        for each in Gg do
            if divide(x_beta, Groebner:-LeadingMonomial(each, ordering_)) then 
            # if x_beta is divisible by any leading monomial in current GB
            # add the S-polynomial component constructed from x_beta and this poly (each) from GB
                Hh := [op(Hh), each * x_beta/Groebner:-LeadingMonomial(each, ordering_)]:
                break:
            end if:
        end do:
        # printf(`\t\t\t\tprocessing time: %a\n`, time()-start):
        start := time():
        Mon_ := GetMonomials(Hh, ordering_): # update list of monomials
    od:
    
    # construct a matrix of coefficients of monomials in Hh
    Mm := Matrix(numelems(Hh), numelems(Mon_)):

    for r from 1 to numelems(Hh) do
        for c from 1 to numelems(Mon_) do
            broken_poly := [op(expand(Hh[r]))]:
            for each in broken_poly do
                if type(each/Mon_[c], numeric) then
                    Mm[r,c] := each/Mon_[c]:
                    break:
                else
                    Mm[r,c] := 0:
                end if:
            end do:
        end do:
    end do:
    
    return Mm, Mon_, Hh:
end proc:

Selection := proc(Bb, Ff, ordering_)
    # implements normal selection strategy
    local degrees, dd, p:
    degrees := [
        seq(
            degree(
                lcm(Groebner:-LeadingMonomial(Ff[p[1]], ordering_), Groebner:-LeadingMonomial(Ff[p[2]], ordering_))
                ),
        p in Bb)
    ]:
    dd := min(degrees):
    return dd, select(x->degree(lcm(Groebner:-LeadingMonomial(Ff[x[1]], ordering_), Groebner:-LeadingMonomial(Ff[x[2]], ordering_)))=dd, Bb):
end proc:

# GetSelectionTable := proc(Ff, ordering_)
#     local degrees, dd, p:
#     s := numelems(Ff): # s in the Cox, Little, O'Shea book
#     t := s:
#     Bb := []:
#     for j from 1 to s do
#         for i from 1 to j do
#             if i<j then
#                 Bb := [op(Bb), [i,j]]:
#             end if:
#         end do:
#     end do:
#     degrees := [
#         seq(
#             degree(
#                 lcm(Groebner:-LeadingMonomial(Ff[p[1]], ordering_), Groebner:-LeadingMonomial(Ff[p[2]], ordering_))
#                 ),
#         p in Bb)
#     ]:
#     degree_to_pair := table([]):
#     for dd in degrees do
#         lcms := []:
#         Bb_prime := select(x->degree(lcm(Groebner:-LeadingMonomial(Ff[x[1]], ordering_), Groebner:-LeadingMonomial(Ff[x[2]], ordering_)))=dd, Bb):
#         for idxs in Bb_prime do
#             i := idxs[1]:
#             j := idxs[2]:
#             LM_i := Groebner:-LeadingMonomial(Ff[i], ordering_):
#             LM_j := Groebner:-LeadingMonomial(Ff[j], ordering_):
#             LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Ff[i], ordering_)[2]*Ff[i]:
#             lcms := [op(lcms), LCM]:
#             LCM := lcm(LM_i, LM_j)/Groebner:-LeadingTerm(Ff[j], ordering_)[2]*Ff[j]:
#             lcms := [op(lcms), LCM]:
#         end do:
#         # printf("time for selection: %a\n", time()-start):
#         lcms := [op({op(lcms)})]:
#         if assigned(degree_to_pair[dd]) then
#             degree_to_pair[dd] := degree_to_pair[dd] union {seq(Groebner:-LeadingMonomial(each, ordering_), each in lcms)}:
#         else
#             degree_to_pair[dd] := {seq(Groebner:-LeadingMonomial(each, ordering_), each in lcms)}:
#         end if:
#     end do:
#     return degree_to_pair:
# end proc:
