
dofile('hcell.lua');

-- default
Lambda = 15 -- mM

V0    = -40 -- mV


gamma = 4.0; -- mV/mM
kQ    = 0.05;
Qout  = 10;
Qini  = 140;
sigma  = 1.00377242943013

Aini  = 20
Aout  = 140



gamma    = 3.64007025527833
k7       = 0.80675570443902
kQ       = 0.0489654980523191
kA_ratio = 0.1

sigma  = 1.00377242943013

passive = "passive.lua"
phi     = 1

if y_file_exists(passive) then
print("Loading " .. passive .. "..." )
dofile(passive)
print("...done")
end
