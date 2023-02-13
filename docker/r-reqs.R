options(
repos = c(
  ggseg = 'https://ggseg.r-universe.dev',
  CRAN = 'https://cloud.r-project.org') 
#lib = '/opt/conda/lib/R/library'
)

package_list <- c(
	'sf',
	'ggseg',
	# 'ggsegExtra',
	'ggsegGlasser',
	'ggsegSchaefer',
	'ggsegDesterieux',
	'ggseg3d', 
	'ggrepel',
	'ggridges',
        'ggwordcloud',
	'ggtext',
	'ggpubr',
	'ggh4x',
	'lemon',
	'patchwork',
	'png',
	'pals', 
	'scales',  
	'shades',
	'magick', 
	'png',
	'gdtools',
	'flextable',
	'vctrs=0.5.0'
)
install.packages(package_list)

devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
