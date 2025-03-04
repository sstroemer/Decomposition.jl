include("cache.jl")

Base.Broadcast.broadcastable(model::Benders.DecomposedModel) = Ref(model)

# Dualization is just flagged, and never directly applied as separate step.
apply!(model::Benders.DecomposedModel, attribute::Solver.DualizeModel) = true

# ... same for these.
apply!(model::Benders.DecomposedModel, attribute::Benders.Config.ModelDebug) = true
apply!(model::Benders.DecomposedModel, attribute::Benders.Config.ModelDirectMode) = true

function apply!(model::Benders.DecomposedModel, attribute::Solver.AbstractAttribute)
    models = attribute.model == :main ? [Benders.main(model)] : Benders.subs(model)
    
    ret = true
    for m in models
        ret &= Solver._modify_jump(m, attribute)
    end
    
    return ret
end

function finalize!(model::Benders.DecomposedModel)
    verbosity = get_attribute(model, Benders.Config.ModelVerbosity).verbosity

    if verbosity >= 4
        Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Debug))
    elseif verbosity >= 3
        Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Info))
    elseif verbosity >= 1
        Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Warn))
    else
        Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Error))
    end

    # Construct & apply all annotations if this is called for the first time.
    if isempty(model.models)
        @timeit model.timer "annotations (generate)" annotate!(model)
        @timeit model.timer "annotations (apply)" apply_annotations!(model)
    end

    # Apply all attributes, that are flagged as dirty.
    for ac in model.attributes
        ac.dirty || continue

        # Skip all "not-to-apply" attributes.
        (ac.attribute isa Benders.Config.ModelVerbosity) && continue

        ret = @timeit model.timer "apply attributes" apply!(model, ac.attribute)
        if !ret
            @error "Failed to apply attribute" model_type = typeof(model) attribute = ac.attribute
        else
            ac.dirty = false
        end
    end

    return nothing
end
