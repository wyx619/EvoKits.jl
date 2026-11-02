# Labels belong to CommunityData; overriding them here could silently associate
# abundance columns with different species or results with different sites.
function _community_data_keywords(kwargs)
    for key in (:tip_labels, :sample_labels)
        haskey(kwargs, key) && throw(ArgumentError("CommunityData owns its labels; set $key when constructing the data"))
    end
    return nothing
end

function _community_tree_data(data, tree)
    if data.species_labels == tree.tip_labels
        _validate_community_data(data)
        return data
    end
    return align_community(data, tree)
end

function weighted_endemism(data::CommunityData; kwargs...)
    _community_data_keywords(kwargs)
    _validate_community_data(data)
    return weighted_endemism(data.matrix; sample_labels = data.sample_labels, kwargs...)
end

function beta_diss(data::CommunityData; kwargs...)
    _community_data_keywords(kwargs)
    _validate_community_data(data)
    return beta_diss(data.matrix; sample_labels = data.sample_labels, kwargs...)
end

# These existing result types carry site labels. Alignment retains the full
# tree, supplies zero columns for unrepresented tips, and owns its arrays.
for metric in (:pd, :pd_dataframe, :phylo_endemism, :relative_pd, :phylobeta,
    :ses_pd, :ses_pd_dataframe, :ses_mpd, :ses_mpd_dataframe,
    :ses_mntd, :ses_mntd_dataframe, :prepare_community)
    @eval function $metric(data::CommunityData, tree::CompactTree; kwargs...)
        _community_data_keywords(kwargs)
        aligned = _community_tree_data(data, tree)
        return $metric(aligned.matrix, tree; sample_labels = aligned.sample_labels, kwargs...)
    end
end

# Matrix/vector results retain the input row order; their existing result
# types are unchanged. tip_weights in pd_bootstrap keep tree-tip order.
for metric in (:unifrac, :phylosor, :comdist, :comdistnt, :raoD,
    :psv, :psr, :pse, :psc, :psd, :pd_calc, :pd_bootstrap)
    @eval function $metric(data::CommunityData, tree::CompactTree; kwargs...)
        _community_data_keywords(kwargs)
        aligned = _community_tree_data(data, tree)
        return $metric(aligned.matrix, tree; kwargs...)
    end
end
