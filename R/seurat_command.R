# script of slightly modified versions of the seurat command script in Seurat
# commands should be added to the seurat@command attribute if the seurat object
# is modified. 

# function needed for make_screp_seurat_cmd
#' @keywords internal
seurat_extractfield <- function(string, field = 1, delim = "_") {
	fields <- as.numeric(
		x = unlist(x = strsplit(x = as.character(x = field), split = ","))
	)
	if (length(x = fields) == 1) {
		return(strsplit(x = string, split = delim)[[1]][field])
	}
	return(paste(
		strsplit(x = string, split = delim)[[1]][fields], collapse = delim
	))
}

# seurat's command adding but if a param is a dataframe or list of dataframes,
# completely omits them to save memory. 
#' @keywords internal
make_screp_seurat_cmd <- function(call_time, assay) {
	
	if (as.character(x = sys.calls()[[1]])[1] == "do.call") {
		call_string <- deparse(expr = sys.calls()[[1]])
		command_name <- as.character(x = sys.calls()[[1]])[2]
	} else {
		command_name <- as.character(
			x = deparse(expr = sys.calls()[[sys.nframe() - 1]])
		)
		command_name <- gsub(
			pattern = "\\.Seurat",
			replacement = "",
			x = command_name
		)
		call_string <- command_name
		command_name <- seurat_extractfield(
			string = command_name,
			field = 1,
			delim = "\\("
		)
	}
	
	argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
	argnames <- grep(
		pattern = "object",
		x = argnames,
		invert = TRUE,
		value = TRUE
	)
	argnames <- grep(
		pattern = "anchorset",
		x = argnames,
		invert = TRUE,
		value = TRUE
	)
	argnames <- grep(
		pattern = "\\.\\.\\.",
		x = argnames,
		invert = TRUE,
		value = TRUE
	)
	
	params <- list()
	p.env <- parent.frame(n = 1)
	argnames <- intersect(x = argnames, y = ls(name = p.env))
	for (arg in argnames) {
		param_value <- get(x = arg, envir = p.env)
		if (is_seurat_object(param_value) || is_df_or_list_of_df(param_value)) {
			next
		}
		params[[arg]] <- param_value
	}
	
	command_name <- sub(
		pattern = "[\\.]+$",
		replacement = "",
		x = command_name,
		perl = TRUE
	)
	command_name <- sub(
		pattern = "\\.\\.", replacement = "\\.", x = command_name, perl = TRUE
	)
	
	# return the command object
	methods::new(
		Class = 'SeuratCommand',
		name = command_name,
		params = params,
		time.stamp = call_time,
		call.string = call_string,
		assay.used = assay
	)
}
