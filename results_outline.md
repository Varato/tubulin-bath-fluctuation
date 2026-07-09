# Paper writing outline


### [Metadata]

### Paper title candidates:

1. Multi-timescale noise forbids long-range coherent exciton transport in tubulin tryptophan networks
2. Fast solvent fluctuations enable incoherent hopping but not coherent exciton delocalization in tubulin

### Writing style
- Language: concise. compact. accurate.
- Convey core information, data with tables and figures. Tables and figures must also be compact. If one thing can be shown in one figure, don't make two.
- Terms, symbols must be consistent. Abbreviations can only be used after its first full term appearing.


### Tryptophan namings

```Python
TRP_NAMING = {
    ("A", 21):  ("Trp1", "αW21"),
    ("A", 346): ("Trp2", "αW346"),
    ("A", 388): ("Trp3", "αW388"),
    ("A", 407): ("Trp4", "αW407"),
    ("B", 21):  ("Trp5", "βW21"),
    ("B", 101): ("Trp6", "βW101"),
    ("B", 344): ("Trp7", "βW344"),
    ("B", 397): ("Trp8", "βW397"),
}
```

### Color convention

Fluctuation sources
```Python
SOURCE_COLORS = {
    "protein":    "#7E57C2",
    "water":      "#29B6F6",
    "nucleotide": "#FF7043",
    "ions":       "#f7ef53",
}
```

Tryptophan colors: use when we need colors to distinguish Trps
```Python
# Colour palette for the 8 Trp sites (consistent across figures)
TRP_CMAP = plt.cm.tab10
TRP_COLORS = [TRP_CMAP(i) for i in range(8)]
TRP_COLORS[6] = TRP_CMAP(9)
```



## Results section outline

For data efficiency, two MD trajectories are prepared to capture the site-energy fluctuations at different timescales:

[Use a small table, or describe it in text]
- 40ns @ 10ps: [mention number of frames, Nyquest freq], named slow trajecory later
- 2ns @ 10fm: [mention number of frames, Nyquest freq], named fast trajectory later

The two trajectories have exactly same MD configurations, taken after deep thermal equilibrium, and the underlying integral time step is 2fs. The only difference is the time interval we save data. The site-energy fluctuations are computed based on Linear Stark Effect. The details for the MD setup and compuation on the fluctuaion are articulated in the Method section. Here we focus on the results.

### Magnitude and statistical properties of site-energy fluctuations

Table 1
[A comprehensive table for all key data for later refering. For each Trps, we list its sigma_total, f1, T1, f2, T2, f3 T3, SASA, theta_local, d_nucl. Also a line for mean values over the 8 Trps]

The computed site-energy fluctuations are well normally distributed. The standard deviations extracted from the slow trajectory are denoted $\sigma_total$, and the deviations extracted fro mthe fast trajectory are denoted by $\sigma_fast$. Both are shown in Table 1.

Since the fast trajecory only lasts 2ns, it undersamples slow mode at nanosecond timescale, $sigma_fast$ is on average smaller then $\sigma_{total}$. We consider $sigma_total$ contains fluctuation mode on all timescales (hence the notation). $\sigma_{total}/J$ ranges from [] to [], indicating high discorder level that should supress coherent energy transfer (Anderson Localization) within the network (discussed later).


### Site-energy fluctuations exhibit three well-separated relaxation timescales

Fig 1
[Stiched autocorrelations and PSD]

Autocorrelations of the site-energy fluctuations are computed from both the fast and slow MD trajectories, denoted by $A_m^{\text{fast}} and A_m^{\text{slow}}$. By stiching them, a multiresolution ACF and Power spectrum are formed (Fig 1).

All ACF arrays in the plot are self-normalized, making the first value (at lag 0) $A_m^{fast}[0] = A_m^{slow}[0] = 1$. The second values, however, correpond to different lag times due to the different sampling rate: A_m^{\text{slow}}[1] is for the lag $Delta t = 10ps$. The stiched ACF arrays exhibit a gap at 10ps because the fast trajectory undersamples slow fluctuation mode, hence has lower deviations during the short lag.

We find the stiched ACFs are best fitted (Fig 1 (a), solid line) by a tri-exponential form
$$
A_m(\tau) = f_1 \mathrm{e}^{-\tau/\T_1} + f_2 \mathrm{e}^{-\tau/\T_2} + f_3 \mathrm{e}^{-\tau/\T_3}
$$
where $f1, f2, f3, T1, T2, T3$ are fitting parameters, shown in Table 1. [Report fitting performence here]. Therefore, the site-energy fluctuations in the tubulin system can be approximated by the superposition of three independent OU processes, with the averaged timescales:

[reporting mean values here, and ]
- T1: 0.044 fs: water libration, f1 = 0.5
- T2: 1.7 ps: hydrogen bond ???, f2 = 0.33
- T3: 2,663 ps: protein structure relaxition, f3 = 0.17
Those values matches the magnitude reported in the literature (Refs needed)


### Source decomposition for site-energy fluctuations and screening effect

Because the site-energy fluctuations are computed using full-atom electric field exerted on the indole ring of each tryptophan, we can investigate the contribution from different sources: water, ions, nucleotides (GTP/GDP attached to the tubulin), and protein.

Fig 2
[left: decomposed sigma, bar chart, shoing sigma_water, sigma_ions, etc. And also sigma_total. And also a direct sum sqrt(sigma^2_g)

right: decomposed PSD for each sources, averaged over all proteins
]

We find that: water contributes most. [put our analysis with numbers from PHASE3]
Also, the fact that sigma_total is [x] times smaller then $sqrt{\sum sigma^2_g}$ indicating the screening effect: the fluctuations cancels each other.

Further, we did the tri-exponential fitting on each source's ACF seperately, and find the matched Ts. This is as a cross-validation


### Impacts of static and dynamic noise on exciton dynamics

Finally, we explore how the above described site-energy fluctuations affects the excitonic dynamics in the tryptophan network. In Craddock 2014, Hanke-Strobl is used to describe a pure dephasing process, which requires Markovian white-noise approximation. Our noise model, however, exhibites color noise at three time scales. Therefore, Monte Carlo approach is applied to directly solve the Stochastic Schrodinger Equation (SSE).

The fitted ACF parameters are used to sample the superposed OU processes for each tryptophan. The fluctuations are applied according to the OU processes to the time-dependent Hamiltonian at a time interval dt = [requires numbers], then a unitary evolutation is performed in the small time windown. $N_{MC} = []$ Monte Carlo trajectories are taken, and the ensemble averaged population dynamics is shown in Fig 3, where the initial state is taken as Trp4 excited.

Fig 3
[Population dynamics]

Interestingly, with only the slowest noise mode $T3$, we see coherent oscilations, but the exciton is strongly localized [numbers needed]. If we add up fast mode at timescales $T1$ and/or $T2$, the oscilation disappears, but the energy can be transfered to the strongly coupled parterner. Notice this transfer is incoherent, and only a ensemble averaged result. This phenomenon is called Environment assisted quantum transport (reference: Environment-assisted quantum transport, 2008). [close sentence needed here]



## Appendix

### Spatial correlations

As a side results, we computed the spatial correlations.

Fig 
[The spatial correlation matrix]

Conclusion is that the site-site spatial correlations are negligible.


### Site-energy fluctuation traces: slow and fast

fig
[4x2 fluctuation plot, with slow + fast. We can show the last 2ns of the slow and the 2ns fast to make both visible and comparable]


### OU process sampling

Equation
$$
x(t + dt) = x(t) * \mathrm{e}^{-dt/T} + \sqrt{1 - \mathrm{e}^{-2dt/T}}\cdot \sigma \xi (t) 
$$
where $T$ and $\sigma^2$ are the corellation time and the variance of the OU process, $\xi(t) \~ N(0, 1)$

Fig 
[the OU sum ACF vs fitted ]

The figure shows that in the long rung, the sum of the three OU process sampleings has identical ACF to our fitted result.

