
function create_tapes(circV, all_strings, W, max_sins, min_abs_coeff, all_treesU, params, nq)
    all_tapes = Vector{Any}(undef, size(all_treesU)[1])
    for j=1:size(all_treesU)[1]
        sigma, site = all_strings[j][2], all_strings[j][1]
        all_tapes[j] = ReverseDiff.compile(ReverseDiff.GradientTape(p->apply_string_to_hybcirc_and_evaluate(all_treesU[j], circV, sigma, site, -1. * p, W, max_sins, min_abs_coeff, nq), [params..., 1.]))
    end 
    return all_tapes
end


function apply_string_to_hybcirc_and_evaluate(treeU, circ, sigma, site, thetas, maxW, maxsins, min_abs_coeff, nq)
    CoeffType = eltype(thetas)
    symb = PauliSum(nq, CoeffType)
    add!(symb, sigma, site, thetas[end])
    wrapped_symb = wrapcoefficients(symb, PauliFreqTracker)
    wrapped_symb = propagate!(circ, wrapped_symb, thetas[1:end-1]; max_weight=maxW, max_freq=Inf, max_sins=maxsins, min_abs_coeff=min_abs_coeff)
    loc_sum_res = overlapwithpaulisum(wrapped_symb, treeU)
    return loc_sum_res
end 

function Cloc(all_strings, params, nq, tapes; return_grad=false)
    all_sum_res = zeros(length(all_strings))
    all_grad_res = zeros(size(params)[1], length(all_strings))

    for index = 1:length(all_strings)
        loc_grad = zeros(size(params)[1]+1)
        ReverseDiff.gradient!(loc_grad, tapes[index], [params..., -1.])
        all_grad_res[:, index] = loc_grad[1:end-1]
        all_sum_res[index] = -loc_grad[end] #the minus comes from the fact that the tape is compiled as p -> -1*p 
    end
    
    if return_grad
        cloc_val = 1 / 2 - (sum(all_sum_res)/(6*nq))
        return abs(cloc_val), sign(cloc_val) * (-1) .* sum(all_grad_res, dims=2) ./ (6*nq)
    end
    return abs(1 / 2 - (sum(all_sum_res)/(6*nq)))
end 