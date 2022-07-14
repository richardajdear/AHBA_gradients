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
	'ggsegExtra',
	'ggseg3d', 
	'ggrepel',
	'ggridges',
	'patchwork',
	'pals', 
	'scales', 
	'ggpubr', 
	'magick', 
	'png',
	'gdtools',
	'flextable'
)
install.packages(package_list)

#install.packages("languageserver", repos='https://cloud.r-project.org/')

devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)
