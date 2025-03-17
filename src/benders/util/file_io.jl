function write_to_file(model::Benders.DecomposedModel; filename::String = "")
    sio = IOBuffer()
    JSON3.pretty(
        sio,
        OrderedDict(
            "meta" => OrderedDict(
                "version" => "0.1.0",
                "algorithm" => "benders",
                "name" => model.name,
                "created" => model.info[:stats][:created],
                "started" => model.info[:stats][:started],
            ),
            "inputs" => OrderedDict(
                "timesteps" => get_attribute(model, Config.TotalTimesteps, :T, nothing),
                "splits" => get_attribute(model, Config.NumberOfTemporalBlocks, :n),
            ),
            "attributes" => showtostr.(model.attributes),
            "models" =>
                OrderedDict("main" => _showtostr(Benders.main(model)), "subs" => _showtostr.(Benders.subs(model))),
            "cuts" => OrderedDict(
                "feasibility" => length(model.cuts[:feasibility]),
                "optimality" => length(model.cuts[:optimality]),
            ),
            "history" => filter.(k -> (k.first != :added_cuts_con), model.info[:history]),
            "log" => join(model.log, "\n"),
        ),
        JSON3.AlignmentContext(; alignment = :Left, indent = 2);
        allow_inf = true,
    )

    json_str = String(take!(sio))
    short_hash = bytes2hex(SHA.sha1(json_str))[1:7]
    filename = isempty(filename) ? "$(model.name)_$(short_hash).djl.json" : filename

    open(normpath(mkpath("out"), filename), "w") do f
        write(f, json_str)
        @info "Saved model" filename
    end

    return nothing
end
