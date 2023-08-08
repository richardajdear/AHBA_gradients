options(
repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org') 
#lib = '/opt/conda/lib/R/library'
)

package_list <- c(
	'httpgd',
	#'sf', # try mamba installing instead
	'Rcpp',
	'ggseg',
	# 'ggsegExtra',
	'ggsegGlasser',
	'ggsegSchaefer',
	'ggsegDesterieux',
	# 'ggseg3d',
	'ggrepel',
	# 'ggridges',
    # 'ggwordcloud',
	'ggtext',
	# 'ggpubr',
	'ggpmisc',
	'gggrid',
	# 'vctrs=0.5.0',
	'ggh4x',
	'lemon',
	'patchwork',
	'png',
	'pals', 
	'scales',
	'shades',
	# 'magick', 
	'png',
	# 'gdtools',
	# 'flextable',
	'eulerr'
)
install.packages(package_list)
