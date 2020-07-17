
using RCall

R"expand.grid.matrix = function(...) data.matrix(expand.grid(...));
    library(plotly)"
xpndgrid = R"expand.grid.matrix"
plotR = R"plot"
linesR = R"lines"
pointsR = R"points"
polygonR = R"polygon"
legendR = R"legend"

function plotly3dsurf(x, y, zmat)
    @rput x y zmat
    R" library('plotly');
    options(browser = 'true');
    plot_ly(type='surface', z = ~zmat, x = ~x, y = ~y);
    "
end

function plotly3dsurf_mh(x, y, zmat, title)
    @rput x y zmat title
    R" library('plotly');
    options(browser = 'false');

    pl = plot_ly(type='surface', z = ~zmat, x = ~x, y = ~y) %>%
    layout(title=title);

    Sys.setenv('plotly_username'='mheiner');
    Sys.setenv('plotly_api_key'='nyP5pZpV9egI4KcmZsAj');
    api_create(pl, filename=title);
    "
end

function plotly3dsurf_bands(x, y, zmat, zmat2, zmat3, opacity_bands)
    @rput x y zmat zmat2 zmat3 opacity_bands
    R" library('plotly');
    options(browser = 'true');
    plot_ly(type='surface', z = ~zmat, x = ~x, y = ~y, showscale=FALSE) %>%
    add_surface(z = ~zmat2, opacity=opacity_bands) %>%
    add_surface(z = ~zmat3, opacity=opacity_bands);
    "
end


function plotly3dsurf_bands_mh(x, y, zmat, zmat2, zmat3, opacity_bands, title)
    @rput x y zmat zmat2 zmat3 opacity_bands title
    R" library('plotly');
    options(browser = 'false');

    pl = plot_ly(type='surface', z = ~zmat, x = ~x, y = ~y, showscale=FALSE) %>%
    add_surface(z = ~zmat2, opacity=opacity_bands) %>%
    add_surface(z = ~zmat3, opacity=opacity_bands);

    Sys.setenv('plotly_username'='mheiner')
    Sys.setenv('plotly_api_key'='nyP5pZpV9egI4KcmZsAj')
    api_create(pl, filename=title );
    "
end


function plotly3dscat(x, y, z)
    @rput x y z
    R" library('plotly');
    options(browser = 'true');
    plot_ly(z = ~z, x = ~x, y = ~y)"
end

function plotly3dscat_mh(x, y, z, title)
    @rput x y z title
    R" library('plotly');
    options(browser = 'false');
    pl = plot_ly(type='scatter3d', z = ~z, x = ~x, y = ~y, mode='markers') %>%
    layout(title=title);

    Sys.setenv('plotly_username'='mheiner');
    Sys.setenv('plotly_api_key'='nyP5pZpV9egI4KcmZsAj');
    api_create(pl, filename=title );
    "
end
