options(
repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org') 
#lib = '/opt/conda/lib/R/library'
)

package_list <- c(
	#'vctrs',
	#'ellipsis',
	'ggseg',
	# 'ggsegExtra',
	'ggsegGlasser',
	'ggsegSchaefer',
	'ggsegDesterieux',
	'ggseg3d', 
	'ggrepel',
	'ggridges',
	'ggtext',
	'ggpubr',
	'lemon',
	'patchwork',
	'pals', 
	'scales',  
	'shades',
	'magick', 
	'png',
	'gdtools',
	'flextable'
)
install.packages(package_list)

devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
