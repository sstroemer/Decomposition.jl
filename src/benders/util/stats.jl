best_upper_bound(model::Benders.DecomposedModel) = minimum(it[:upper_bound] for it in model.info[:history]; init=+Inf)
best_lower_bound(model::Benders.DecomposedModel) = maximum(it[:lower_bound] for it in model.info[:history]; init=-Inf)

best_gap_abs(model::Benders.DecomposedModel) = gap_abs(best_lower_bound(model), best_upper_bound(model))
best_gap_rel(model::Benders.DecomposedModel) = gap_rel(best_lower_bound(model), best_upper_bound(model))
