module Solution

export SwitchingSolutionView

using ..Types

# ============================================================
# WRAPPER
# ============================================================

struct SwitchingSolutionView{S}
    sol::S
end


# ============================================================
# INTERNAL ACCESSORS
# ============================================================

tau(s::SwitchingSolutionView)    = s.sol.tau
before(s::SwitchingSolutionView) = s.sol.solution_before
after(s::SwitchingSolutionView)  = s.sol.solution_after
switched(s::SwitchingSolutionView) = s.sol.switched


# ============================================================
# PROPERTY ACCESS (PUBLIC API)
# ============================================================

function Base.getproperty(s::SwitchingSolutionView, name::Symbol)

    # --- exposed fields ---
    name === :tau      && return tau(s)
    name === :before   && return before(s)
    name === :after    && return after(s)
    name === :switched && return switched(s)

    # --- allow access to underlying struct if needed ---
    name === :sol && return getfield(s, :sol)

    return getfield(s, name)
end


# ============================================================
# TIME QUERY INTERFACE
# ============================================================

function (s::SwitchingSolutionView)(t)
    sol = s.sol

    if !sol.switched || t < sol.tau
        return sol.solution_before(t; idxs=1)
    else
        return sol.solution_after(t; idxs=1)
    end
end

end