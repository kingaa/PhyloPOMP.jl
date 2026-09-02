# SI2R (superspreading) model declaration.
#
# Compartments: S (susceptible), I_L (low-rate spreader), I_H (high-rate
# super-spreader), R (recovered).
# Demes: I_L, I_H  (two infectious compartments, matching the SI2R QMD's
#   D = {L, H}).
#
# Jump marks (9 total):
#   TL  — transmission by low-rate spreader, BIRTH, r=(2,0), from=I_L
#   TH  — transmission by super-spreader, BIRTH, r=(1,1), from=I_H
#   L   — transition to super-spreading, MIGRATION I_L→I_H, r=(0,1)
#   H   — transition to low-rate, MIGRATION I_H→I_L, r=(1,0)
#   RL  — recovery of low-rate spreader, DEATH from I_L
#   RH  — recovery of super-spreader, DEATH from I_H
#   W   — waning of immunity, NEUTRAL
#   SL  — sampling of low-rate spreader, SAMPLE from I_L
#   SH  — sampling of super-spreader, SAMPLE from I_H
#
# Reference: si2r_model.qmd (independently hand-derived filter equation).

@mgp SI2R begin
    compartments = (S, I_L, I_H, R)
    demes = (I_L, I_H)
    params = (β, κ, γ, ω, ψ, η_L, η_H, N)

    # TL: transmission by a low-rate spreader. Parent in I_L, both offspring
    # in I_L (same-deme fork). r = (2,0).
    @event TL rate=β*S*I_L/N  pop=(S=-1, I_L=+1) move=fork(I_L => I_L, I_L) kind=regular

    # TH: transmission by a super-spreader. Parent in I_H, offspring go to
    # I_L (the new exposed enters low-rate) and I_H (parent stays). r = (1,1).
    @event TH rate=κ*β*S*I_H/N  pop=(S=-1, I_L=+1) move=fork(I_H => I_L, I_H) kind=regular

    # L: transition to super-spreading behavior. MIGRATION I_L → I_H. r = (0,1).
    @event L rate=η_L*I_L  pop=(I_L=-1, I_H=+1) move=swap(I_L => I_H) kind=regular

    # H: transition to low-rate spreading. MIGRATION I_H → I_L. r = (1,0).
    @event H rate=η_H*I_H  pop=(I_H=-1, I_L=+1) move=swap(I_H => I_L) kind=regular

    # RL: recovery of low-rate spreader. DEATH from I_L.
    @event RL rate=γ*I_L  pop=(I_L=-1, R=+1) move=chop(I_L) kind=regular

    # RH: recovery of super-spreader. DEATH from I_H.
    @event RH rate=γ*I_H  pop=(I_H=-1, R=+1) move=chop(I_H) kind=regular

    # W: waning of immunity. NEUTRAL (no effect on genealogy).
    @event W rate=ω*R  pop=(R=-1, S=+1) move=none kind=regular

    # SL: sampling of low-rate spreader. SAMPLE from I_L.
    @event SL rate=ψ*I_L  pop=() move=sample(I_L) kind=singular

    # SH: sampling of super-spreader. SAMPLE from I_H.
    @event SH rate=ψ*I_H  pop=() move=sample(I_H) kind=singular
end

# we are assuming the integer vector is added to the population vector since the state can be more complicated than a vector of integers. 