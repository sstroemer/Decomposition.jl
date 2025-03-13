#!/usr/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_02.jl
julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_03.jl
julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_05.jl
julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_06.jl
julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_07.jl
julia --project=$SCRIPT_DIR $SCRIPT_DIR/analyse_08.jl
