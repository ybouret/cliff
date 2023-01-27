dofile('pcell.lua')


E0    = 30
K5    = 20
f7    = 30
Z5    = 0.0

phi    = 1
decay  = 9
eta_f  = 1.0136
eta_r  = 1.0095


active = "active.lua"

if y_file_exists(active) then
print("Loading " .. active .. "..." )
-- dofile(active)
print("...done")
end

delta7 =
{
    -- "K5", "f7", "K5:f7", "eta_f", "f7:eta_f", "eta_r", "f7:eta_r", "f7:eta_f:eta_r", "K5:f7:eta_f:eta_r"

  "f7", "eta_f", "eta_f:f7", "eta_r", "eta_r:f7" , "eta_r:f7:eta_f"
}

intake =
{
    "phi"
}

loops=1

-- set y2tics;
-- set logscale x;
-- plot [1:3600] [0:15] '../../data/nhe1_delta7_15mM_v3.txt' w lp, 'nhe1_delta7_15mM_v3.fcn.txt' u 1:7 w l, '../../data/nhe1_intake_15mM_v3.txt' w lp axis x1y2, 'nhe1_intake_15mM_v3.fcn.txt' u 1:8 w l axis x1y2, 'nhe1_intake_15mM_v3.fit.txt' u 1:3 w lp axis x1y2
