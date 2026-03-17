
using PhyloPOMP

@info "Newick format reading"

x1 = readlines("MERS_274_sCoal_phylopomp.nwk")

@time parse_newick(x1, 0.0)

x2 = [
	"(((((b_1_1:0.1,s_0_2:0.1)g_0_1:0.2,b_0_3:0.3)g_0_2:0.4,b_1_4:0.5)g_0_3:0.6)b_0_8:0.2,b_0_5:0.7)m_0_4:0.8;",
	"(b_0_6:0.1,b_0_7:0.1)m_0_5:0.5;",
]

@time parse_newick(x2, 5.0)
