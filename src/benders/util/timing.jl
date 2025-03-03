total_wall_time(model::Benders.DecomposedModel) = sum(it[:time][:wall] for it in model.info[:history]; init=0)
total_cpu_time(model::Benders.DecomposedModel) = sum(it[:time][:cpu] for it in model.info[:history]; init=0)

get_main_time(model::Benders.DecomposedModel) = TimerOutputs.time(model.timer["main"])
get_sub_time(model::Benders.DecomposedModel; index::Int) = TimerOutputs.time(model.timer["sub"]["[$(index)]"])

function summarize_timings(model::Benders.DecomposedModel)
    to = TimerOutputs.TimerOutput()
    
    to.inner_timers["main"] = model.timer["main"]

    to.inner_timers["sub"] = TimerOutputs.merge(values(model.timer["sub"].inner_timers)...)
    to.inner_timers["sub"].name = "sub [1 ... $(length(Benders.subs(model)))]"
    
    return to
end
