#two.py    - dash app
import time
import scanpy as sc
import numpy as np
import pandas as pd
import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import plotly.express as px

# Bootstrap theme for styling
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "Single Cell Expression Atlas for Autonomic Nervous System"

# Load two datasets
NEURON_FILE = "data/2_sym_para_master_neurons.h5ad"
ALLCELLS_FILE = "data/2_sym_para_master_allcells.h5ad"

adata_neurons = sc.read_h5ad(NEURON_FILE)
adata_all = sc.read_h5ad(ALLCELLS_FILE)

print("Loaded neuron dataset:", adata_neurons.shape)
print("Loaded all cells dataset:", adata_all.shape)

# Build base dataframes from the AnnData objects
def build_base_dfs(adata):
    umap_orig = adata.obsm["X_umap_no_harmony"]
    base_df_orig = pd.DataFrame({
        "UMAP1": umap_orig[:, 0],
        "UMAP2": umap_orig[:, 1],
        "UMAP3": umap_orig[:, 2],
        "Ganglia": adata.obs["batch2"].values,
        "PMID": adata.obs["PMID"].values,
        "ANS_Division": adata.obs["ANS_Division"].values,
        "Dataset": adata.obs["dataset"].values
    }, index=adata.obs_names)
    
    umap_harm = adata.obsm["X_umap_harmony"]
    base_df_harm = pd.DataFrame({
        "UMAP1": umap_harm[:, 0],
        "UMAP2": umap_harm[:, 1],
        "UMAP3": umap_harm[:, 2],
        "Ganglia": adata.obs["batch2"].values,
        "PMID": adata.obs["PMID"].values,
        "ANS_Division": adata.obs["ANS_Division"].values,
        "Dataset": adata.obs["dataset"].values
    }, index=adata.obs_names)
    return base_df_orig, base_df_harm

# Create dataframes for neurons and all cells
base_df_original_neurons, base_df_harmony_neurons = build_base_dfs(adata_neurons)
base_df_original_all, base_df_harmony_all = build_base_dfs(adata_all)

# Set the static HTML for the app (index_string)
app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        <meta charset="UTF-8">
        <title>Single Cell Expression Atlas for Autonomic Nervous System</title>
        <link href="https://fonts.googleapis.com/css2?family=Open+Sans:wght@400;700&display=swap" rel="stylesheet">
        <style>
            html { scroll-behavior: smooth; }
            body {
                font-family: 'Open Sans', sans-serif;
                background-color: #f9f9f9;
                margin: 0; padding: 0; color: #333;
            }
            h2 {
                color: #333;
                margin-bottom: 20px;
            }
            .top-nav {
                position: absolute; top: 10px; right: 10px; z-index: 1000;
            }
            .top-nav a {
                margin-left: 15px; color: #2196F3; text-decoration: none; font-weight: 600;
            }
            .top-nav a:hover {
                text-decoration: underline;
            }
            .segmented-control {
                display: flex; flex-direction: row; border: 1px solid #ccc;
                border-radius: 5px; overflow: hidden; width: fit-content;
            }
            .segmented-control .form-check-inline {
                display: inline-flex !important; margin-right: 0;
            }
            .segmented-control .form-check-label {
                text-align: center; padding: 10px 20px; cursor: pointer;
                border-right: 1px solid #ccc; margin: 0;
            }
            .segmented-control .form-check-label:last-child {
                border-right: none;
            }
            .segmented-control .form-check-input {
                display: none;
            }
            .segmented-control .form-check-input:checked + label {
                background-color: #2196F3; color: #fff;
            }
            /* Gene search + submit button */
            input#gene-input {
                border: 1px solid #aaa; border-radius: 5px; padding: 8px;
                transition: border-color 0.2s, box-shadow 0.2s;
            }
            input#gene-input:focus {
                border-color: #2196F3; outline: none;
                box-shadow: 0 0 5px rgba(33, 150, 243, 0.5);
            }
            button#submit-gene {
                font-family: 'Open Sans', sans-serif;
                border: 1px solid #aaa; border-radius: 5px;
                padding: 8px 12px; background-color: #fff; cursor: pointer;
                transition: background-color 0.2s, border-color 0.2s, box-shadow 0.2s;
            }
            button#submit-gene:hover {
                background-color: #f0f0f0;
            }
            button#submit-gene:focus {
                border-color: #2196F3; box-shadow: 0 0 5px rgba(33, 150, 243, 0.5);
                outline: none;
            }
            /* Subtitle link */
            .subtitle-link {
                color: #2196F3; text-decoration: underline;
            }
            .subtitle-link:hover {
                color: grey;
            }
            /* Add a border to each individual checkbox item */
            .bordered-checklist .form-check {
                border: 1px solid #ccc;
                border-radius: 5px;
                margin: 5px 5px 5px 0;
                padding: 3px 6px;
                display: inline-block;
            }
            /* Style disabled checkbox labels as gray */
            .bordered-checklist .form-check input:disabled + label {
                color: #aaa;
            }
        </style>
        <script>
            document.addEventListener("DOMContentLoaded", function() {
                document.title = "Single Cell Expression Atlas for Autonomic Nervous System";
                Object.defineProperty(document, "title", {
                    set: function(newTitle) {},
                    get: function() { return "Single Cell Expression Atlas for Autonomic Nervous System"; }
                });
            });
        </script>
        <link rel="icon" href="data:,">
    </head>
    <body>
        <!--This Expression Atlas was built by Avedis Tufenkjian, 2025 -->
        <div class="top-nav">
            <a href="https://www.nature.com/articles/s41583-025-00941-2" target="_blank">Reference</a>
            <a href="#about-section">About</a>
            <a href="https://github.com/YTwTJ/Molecular-and-Functional-Diversity-of-the-Autonomic-Nervous-System" target="_blank">Source Code</a>
            <a href="https://okalab.caltech.edu/" target="_blank">Lab Website</a>
        </div>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""

## Layout and styling for dropdowns
style_filter_box = {
    "display": "inline-block",
    "marginRight": "20px",
    "verticalAlign": "top",
    "border": "1px solid #ccc",
    "padding": "10px",
    "borderRadius": "6px"
}

app.layout = html.Div([
    html.H2("Single Cell Expression Atlas for Autonomic Nervous System",
            style={"textAlign": "center", "marginBottom": "10px"}),

    html.Div([
        html.Span("An interactive platform for exploring publicly available transcriptomic datasets",
            style={"fontSize": "16px", "color": "#333"}),
        
        html.Br(),  
        html.Span(
            [
                "Read our Review for more information ",
                html.A("here", href="https://rdcu.be/euEhn", target="_blank",style={"color": "#2196F3", "fontWeight": "bold"}),
            ],
            style={"fontSize": "16px", "color": "#333"},
        ),
        
    ], style={"textAlign": "center", "marginBottom": "30px"}),

    # Cascading filtering options
    html.Div([
        html.Div([
            html.Label("1) Select ANS Division:", style={"fontWeight": "bold"}),
            dbc.Checklist(
                id='ans-checklist',
                options=[],
                value=[],
                inline=True,
                className="bordered-checklist" 
            )
        ], style=style_filter_box),

        html.Div([
            html.Label("2) Select Ganglia:", style={"fontWeight": "bold"}),
            dbc.Checklist(
                id='ganglia-checklist',
                options=[],
                value=[],
                inline=True,
                className="bordered-checklist"
            )
        ], style=style_filter_box),

        html.Div([
            html.Label("3) Select PMID:", style={"fontWeight": "bold"}),
            dbc.Checklist(
                id='pmid-checklist',
                options=[],
                value=[],
                inline=True,
                className="bordered-checklist"
            )
        ], style=style_filter_box)
    ], style={"marginBottom": "20px", "textAlign": "center"}),

    # Combined row with Cell Type Toggle, Gene Search Bar, and UMAP Toggle below filtering options
    html.Div([
        # Cell type toggle on the left
        html.Div(
            dbc.RadioItems(
                id="data-type-toggle",
                options=[
                    {"label": "Neurons", "value": "neurons"},
                    {"label": "All Cells", "value": "all"}
                ],
                value="neurons",
                inline=True
            ),
            className="segmented-control",
            style={"marginRight": "20px"}
        ),
        # Gene search bar in the center
        html.Div([
            html.Label("Gene Name:"),
            dcc.Input(
                id="gene-input", type="text", placeholder="e.g. 'Th'",
                list="gene-list", style={"width": "200px", "marginLeft": "10px"}
            ),
            html.Datalist(id="gene-list", children=[]),
            html.Button("Submit", id="submit-gene", n_clicks=0,
                        style={"marginLeft": "10px", "backgroundColor": "#fff"})
        ], style={"textAlign": "center", "marginRight": "20px"}),
        # UMAP toggle on the right
        html.Div(
            dbc.RadioItems(
                id="umap-toggle",
                options=[
                    {"label": "Original UMAP", "value": "original"},
                    {"label": "Harmonized UMAP", "value": "harmony"}
                ],
                value="original",
                inline=True
            ),
            className="segmented-control"
        )
    ], style={
        "display": "flex",
        "justifyContent": "center",
        "alignItems": "center",
        "marginBottom": "10px"
    }),

    # Plots and layout
    html.Div([
        html.Div([
            dcc.Graph(
                id="umap-3d-plot",
                style={"height": "900px", "width": "100%"},
                config={"displaylogo": False}
            )
        ], style={"width": "50%", "display": "inline-block", "verticalAlign": "top"}),

        html.Div([
            dcc.Graph(
                id="umap-dataset-plot",
                style={"height": "900px", "width": "100%"},
                config={"displaylogo": False}
            )
        ], style={"width": "50%", "display": "inline-block", "verticalAlign": "top"})
    ], style={"width": "100%", "margin": "auto", "textAlign": "center"}),

    # About section
    html.Div([
        html.H2("About", style={"textAlign": "left", "marginLeft": "20px"}),
        dcc.Markdown("This Expression Atlas of the Autonomic Nervous System provides an interactive platform to visualize cellular heterogeneity and explore gene expression patterns in publicly available transcriptomic datasets.", style={"textAlign": "left", "marginLeft": "20px"})
    ], id="about-section", style={"margin": "20px 0"}),

    html.Div([
        html.H2("Key Features", style={"textAlign": "left", "marginLeft": "20px"}),
        dcc.Markdown("""**Customizable dataset views.** 
                     Three dynamic filtering parameters to display datasets of interest: ANS divisions, specific ganglia or nuclei, and associated publications. The current data includes over 70k cells from eight datasets. 

**Interactive UMAP visualizations.** 
Explore single-cell data through dynamic and zoomable UMAP plots. Compare gene expression profiles and cell dataset origins side-by-side.

**Original and integrated UMAP embeddings.** 
                     Toggle between the original UMAP embedding and Harmony-integrated UMAP embedding. Harmony integration reduces batch-specific biases, aligning datasets to reveal biologically relevant clusters.

Reference our review article: “Molecular and Functional Diversity of the Autonomic Nervous System.” 
""", style={"textAlign": "left", "marginLeft": "20px"})
    ], id="methods-section", style={"margin": "20px 0"}),

    # Dummy stores for camera sync, scroll trigger, and camera reset
    dcc.Store(id="dummy-store-l2r"),
    dcc.Store(id="dummy-store-r2l"),
    dcc.Store(id="scroll-trigger"),
    dcc.Store(id="reset-camera")
], style={"padding": "20px"})

# Callback: Update ANS checklist based on data type toggle
@app.callback(
    [Output("ans-checklist", "options"),
     Output("ans-checklist", "value")],
    [Input("data-type-toggle", "value")]
)
def update_ans_checklist(data_type):
    if data_type == "neurons":
        full_df = base_df_original_neurons
    else:
        full_df = base_df_original_all
    unique_ans = sorted(full_df["ANS_Division"].unique())
    options = [{"label": a, "value": a, "disabled": False} for a in unique_ans]
    return options, unique_ans

# Callback: Update Ganglia checklist based on selected ANS division and data type
@app.callback(
    [Output("ganglia-checklist", "options"),
     Output("ganglia-checklist", "value")],
    [Input("ans-checklist", "value"),
     Input("data-type-toggle", "value")]
)
def update_ganglia_checklist(selected_ans, data_type):
    if data_type == "neurons":
        full_df = base_df_original_neurons
    else:
        full_df = base_df_original_all

    full_ganglia = sorted(full_df["Ganglia"].unique())
    if not selected_ans:
        available_ganglia = []
    else:
        available_ganglia = sorted(full_df[full_df["ANS_Division"].isin(selected_ans)]["Ganglia"].unique())
    
    options = [
        {"label": g, "value": g, "disabled": g not in available_ganglia}
        for g in full_ganglia
    ]
    selected_value = [g for g in full_ganglia if g in available_ganglia]
    return options, selected_value

# Callback: Update PMID checklist based on selected ANS division, ganglia, and data type
@app.callback(
    [Output("pmid-checklist", "options"),
     Output("pmid-checklist", "value")],
    [Input("ans-checklist", "value"),
     Input("ganglia-checklist", "value"),
     Input("data-type-toggle", "value")]
)
def update_pmid_checklist(selected_ans, selected_ganglia, data_type):
    if data_type == "neurons":
        full_df = base_df_original_neurons
    else:
        full_df = base_df_original_all

    full_pm = sorted(full_df["PMID"].unique())
    if not selected_ans:
        available_pm = []
    else:
        available_df = full_df[full_df["ANS_Division"].isin(selected_ans)]
        if selected_ganglia:
            available_df = available_df[available_df["Ganglia"].isin(selected_ganglia)]
        available_pm = sorted(available_df["PMID"].unique())

    options = []
    for p in full_pm:
        parts = p.split(", ")
        if len(parts) == 3:
            new_label = f"{parts[2]}, {parts[0]}, {parts[1]}"
        else:
            new_label = p
        options.append({
            "label": new_label,
            "value": p,
            "disabled": p not in available_pm
        })
    options.sort(key=lambda x: x["label"])
    selected_value = [opt["value"] for opt in options if not opt["disabled"]]
    return options, selected_value

# Callback: Update gene name suggestions
@app.callback(
    Output("gene-list", "children"),
    [Input("gene-input", "value"),
     Input("data-type-toggle", "value")]
)
def update_gene_suggestions(input_text, data_type):
    if data_type == "neurons":
        gene_names = adata_neurons.raw.var_names.tolist() if adata_neurons.raw is not None else []
    else:
        gene_names = adata_all.raw.var_names.tolist() if adata_all.raw is not None else []
    if input_text:
        suggestions = [gene for gene in gene_names if gene.lower().startswith(input_text.lower())]
        suggestions = suggestions[:10]
        return [html.Option(value=gene) for gene in suggestions]
    return []

# Callback: Update 3D UMAP plot based on filtering and gene search
@app.callback(
    Output("umap-3d-plot", "figure"),
    [
        Input("submit-gene", "n_clicks"),
        Input("gene-input", "n_submit"),
        Input("ganglia-checklist", "value"),
        Input("pmid-checklist", "value"),
        Input("ans-checklist", "value"),
        Input("umap-toggle", "value"),
        Input("data-type-toggle", "value")
    ],
    [State("gene-input", "value")]
)
def update_umap(n_clicks, n_submit, selected_ganglia, selected_pmids,
                selected_ans, umap_toggle, data_type, gene_name):
    
    data_type_str = "neurons" if data_type == "neurons" else "all"
    umap_type_str = "Original UMAP" if umap_toggle == "original" else "Harmonized UMAP"
    
    if data_type_str == "neurons":
        base_df = (base_df_original_neurons if umap_toggle == "original"
                   else base_df_harmony_neurons)
    else:
        base_df = (base_df_original_all if umap_toggle == "original"
                   else base_df_harmony_all)

    filtered_df = base_df[
        (base_df["Ganglia"].isin(selected_ganglia)) &
        (base_df["PMID"].isin(selected_pmids)) &
        (base_df["ANS_Division"].isin(selected_ans))
    ]
    
    if not gene_name:
        # Dummy expression column for blank colorbar
        filtered_df = filtered_df.copy()
        filtered_df["Expression"] = 0
        fig = px.scatter_3d(
            filtered_df, x="UMAP1", y="UMAP2", z="UMAP3",
            color="Expression",
            color_continuous_scale=[(0, "#CCCCCC"), (1, "#CCCCCC")],
            range_color=[0, 1],
            title=f"{umap_type_str} (No gene entered)",
            height=900
        )
        fig.update_traces(marker=dict(size=2, opacity=0.8))
        fig.update_layout(
            autosize=True,
            showlegend=False,
            title_x=0.5,
            coloraxis_colorbar=dict(
                orientation="h",
                yanchor="top", y=-0.25,
                xanchor="center", x=0.5,
                len=0.6, thickness=10
            )
        )
        return fig

    if data_type_str == "neurons":
        gene_list = adata_neurons.raw.var_names.tolist()
    else:
        gene_list = adata_all.raw.var_names.tolist()

    matched_gene = next((g for g in gene_list if g.lower() == gene_name.lower()), None)
    
    if matched_gene is None:
        fig = px.scatter_3d(
            filtered_df, x="UMAP1", y="UMAP2", z="UMAP3",
            title=f"Gene {gene_name} not found!", height=900
        )
        fig.update_traces(marker=dict(size=2, opacity=0.8, color="#CCCCCC"))
        fig.update_layout(autosize=True, title_x=0.5)
        return fig

    if data_type_str == "neurons":
        expr = adata_neurons.raw[:, matched_gene].X
    else:
        expr = adata_all.raw[:, matched_gene].X

    if not isinstance(expr, np.ndarray):
        expr = expr.toarray()
    expr = expr.flatten()
    expr_series = pd.Series(expr, index=base_df.index, name="Expression")
    expr_filtered = expr_series.loc[filtered_df.index]

    plot_df = filtered_df.copy()
    plot_df["Expression"] = expr_filtered.values

    fig = px.scatter_3d(
        plot_df, x="UMAP1", y="UMAP2", z="UMAP3",
        color="Expression", color_continuous_scale="Viridis",
        title=f"{umap_type_str}: {matched_gene} Expression",
        hover_data=["Ganglia", "PMID", "ANS_Division", "Expression"],
        height=900
    )
    fig.update_traces(marker=dict(size=2, opacity=0.8))
    fig.update_layout(
        autosize=True,
        scene=dict(xaxis_title="UMAP 1", yaxis_title="UMAP 2", zaxis_title="UMAP 3"),
        coloraxis_colorbar=dict(
            orientation="h", yanchor="top", y=-0.2,
            xanchor="center", x=0.5, len=0.5
        ),
        title_x=0.5
    )
    return fig

# Callback: Update dataset plot (colored by 'Dataset')
@app.callback(
    Output("umap-dataset-plot", "figure"),
    [
        Input("ganglia-checklist", "value"),
        Input("pmid-checklist", "value"),
        Input("ans-checklist", "value"),
        Input("umap-toggle", "value"),
        Input("data-type-toggle", "value")
    ]
)
def update_dataset_plot(selected_ganglia, selected_pmids, selected_ans,
                        umap_toggle, data_type):
    
    data_type_str = "neurons" if data_type == "neurons" else "all"
    if data_type_str == "neurons":
        base_df = (base_df_original_neurons if umap_toggle == "original"
                   else base_df_harmony_neurons)
    else:
        base_df = (base_df_original_all if umap_toggle == "original"
                   else base_df_harmony_all)

    filtered_df = base_df[
        (base_df["Ganglia"].isin(selected_ganglia)) &
        (base_df["PMID"].isin(selected_pmids)) &
        (base_df["ANS_Division"].isin(selected_ans))
    ]
    
    fig = px.scatter_3d(
        filtered_df,
        x="UMAP1", y="UMAP2", z="UMAP3",
        color="Dataset",
        hover_data=["Ganglia", "PMID", "ANS_Division", "Dataset"],
        title="UMAP Colored by Dataset",
        height=900
    )
    fig.update_traces(marker=dict(size=2, opacity=0.8))
    fig.update_layout(
        autosize=True,
        legend=dict(
            orientation="h",
            y=-0.2,
            x=0.5,
            xanchor="center",
            font=dict(size=16),
            itemsizing="constant",
            itemwidth=30
        ),
        title_x=0.5
    )
    return fig

# Clientside callback: Sync camera from left to right plot
app.clientside_callback(
    """
    function(relayoutData) {
        var defaultCam = {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}};
        var cam;
        if (!relayoutData || Object.keys(relayoutData).length === 0 || !relayoutData.hasOwnProperty('scene.camera')) {
            cam = defaultCam;
        } else {
            cam = relayoutData['scene.camera'];
            if (!cam || Object.keys(cam).length === 0) {
                cam = defaultCam;
            }
        }
        if (window.cameraTimeoutL2R) {
            clearTimeout(window.cameraTimeoutL2R);
        }
        window.cameraTimeoutL2R = setTimeout(function() {
            var targetContainer = document.getElementById('umap-dataset-plot');
            if (targetContainer) {
                var targetGraph = targetContainer.getElementsByClassName('js-plotly-plot')[0];
                if (targetGraph) {
                    try {
                        Plotly.relayout(targetGraph, {'scene.camera': cam});
                    } catch (e) { console.error(e); }
                }
            }
        }, 200);
        return window.dash_clientside.no_update;
    }
    """,
    Output("dummy-store-l2r", "data"),
    Input("umap-3d-plot", "relayoutData")
)

# Clientside callback: Sync camera from right to left plot
app.clientside_callback(
    """
    function(relayoutData) {
        var defaultCam = {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}};
        var cam;
        if (!relayoutData || Object.keys(relayoutData).length === 0 || !relayoutData.hasOwnProperty('scene.camera')) {
            cam = defaultCam;
        } else {
            cam = relayoutData['scene.camera'];
            if (!cam || Object.keys(cam).length === 0) {
                cam = defaultCam;
            }
        }
        if (window.cameraTimeoutR2L) {
            clearTimeout(window.cameraTimeoutR2L);
        }
        window.cameraTimeoutR2L = setTimeout(function() {
            var targetContainer = document.getElementById('umap-3d-plot');
            if (targetContainer) {
                var targetGraph = targetContainer.getElementsByClassName('js-plotly-plot')[0];
                if (targetGraph) {
                    try {
                        Plotly.relayout(targetGraph, {'scene.camera': cam});
                    } catch (e) { console.error(e); }
                }
            }
        }, 200);
        return window.dash_clientside.no_update;
    }
    """,
    Output("dummy-store-r2l", "data"),
    Input("umap-dataset-plot", "relayoutData")
)

# Clientside callback: Scroll to the left plot after gene submission
app.clientside_callback(
    """
    function(n_clicks, n_submit) {
        if ((n_clicks && n_clicks > 0) || (n_submit && n_submit > 0)) {
            var plotElement = document.getElementById("umap-3d-plot");
            if (plotElement) {
                plotElement.scrollIntoView({behavior: "smooth"});
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("scroll-trigger", "data"),
    [Input("submit-gene", "n_clicks"),
     Input("gene-input", "n_submit")]
)

# Clientside callback: Reset the right plot's camera orientation after gene submission
app.clientside_callback(
    """
    function(n_clicks, n_submit) {
        if ((n_clicks && n_clicks > 0) || (n_submit && n_submit > 0)) {
            var defaultCam = {"eye": {"x": 1.25, "y": 1.25, "z": 1.25}};
            var targetContainer = document.getElementById('umap-dataset-plot');
            if (targetContainer) {
                var targetGraph = targetContainer.getElementsByClassName('js-plotly-plot')[0];
                if (targetGraph) {
                    Plotly.relayout(targetGraph, {'scene.camera': defaultCam});
                }
            }
        }
        return window.dash_clientside.no_update;
    }
    """,
    Output("reset-camera", "data"),
    [Input("submit-gene", "n_clicks"), Input("gene-input", "n_submit")]
)

server = app.server

if __name__ == "__main__":
    app.run_server(debug=False, port=80)  # Switch to port 80

