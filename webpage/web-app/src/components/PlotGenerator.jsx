import React, { useState, useCallback, useMemo } from 'react';
import Plot from 'react-plotly.js';
import { BarChart3, Loader, Download, RefreshCw } from 'lucide-react';
import './PlotGenerator.css';

// Mock API for demo - replace with your actual plotApi
const plotApi = {
    generatePlotData: async (gridParams, plotType, potParams) => {
        await new Promise(resolve => setTimeout(resolve, 1000));
        const x = Array.from({ length: 100 }, (_, i) => i / 10);
        const y = x.map(val => Math.sin(val) * Math.exp(-val / 10));
        return {
            data: [{ x, y, type: 'scatter', mode: 'lines' }],
            layout: { title: `${plotType} Potential`, paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0.1)' }
        };
    }
};

// Constants moved outside a component
const PLOT_TYPES = [
    { value: 'Harmonic', label: 'Harmonic' },
    { value: 'Polynomial', label: 'Polynomial' },
    { value: 'Barrier', label: 'Barrier' },
    { value: 'MultiBarrier', label: 'Multiple Barrier'},
    { value: 'Gaussian', label: 'Gaussian' },
    { value: 'MultiGaussian', label: 'Multiple Gaussian' },
    { value: 'SuperGaussian', label: 'Super Gaussian' },
    { value: 'Sinusoidal', label: 'Sinusoidal' },
    { value: 'Morse', label: 'Morse' },
    { value: 'Softcore', label: 'Softcore' },
    { value: 'surface_3d', label: '3D Surface Plot' },
    { value: 'energy_levels', label: 'Energy Level Diagram' },
];

const INITIAL_GRID_PARAMS = { rMin: -0.0, rMax: 10.0, nGrid: 100 };
const INITIAL_POT_PARAMS = { D: 100.0, a: 1.5, r0: 2.0 };

const PLOT_CONFIG = {
    responsive: true,
    displayModeBar: true,
    displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d'],
    toImageButtonOptions: {
        format: 'png',
        filename: 'quantum_plot',
        height: 800,
        width: 1200,
        scale: 2,
    },
};

// Reusable input component
const NumberInput = ({ label, value, onChange, step = "0.1" }) => (
    <div className="input-group">
        <label className="input-label">{label}</label>
        <input
            type="number"
            step={step}
            value={value}
            onChange={(e) => onChange(parseFloat(e.target.value))}
            className="input-field"
        />
    </div>
);

// Parameter configurations
const PARAM_CONFIGS = {
    Morse: [
        { key: 'D', label: 'Dissociation Energy (D)', step: '10' },
        { key: 'a', label: 'Width Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Equilibrium Distance (r₀)', step: '0.1' },
    ],
    Softcore: [
        { key: 'D', label: 'Charge (q)', step: '1' },
        { key: 'a', label: 'Softcore Parameter (a)', step: '0.1' },
        { key: 'r0', label: 'Center (r₀)', step: '0.1' },
    ],
};

export default function PlotGenerator() {
    const [gridParams, setGridParams] = useState(INITIAL_GRID_PARAMS);
    const [potParams, setPotParams] = useState(INITIAL_POT_PARAMS);
    const [plotType, setPlotType] = useState('Morse');
    const [isGenerating, setIsGenerating] = useState(false);
    const [plotData, setPlotData] = useState(null);
    const [error, setError] = useState(null);

    // Memoized update functions
    const updateGridParam = useCallback((key, value) => {
        setGridParams(prev => ({ ...prev, [key]: value }));
    }, []);

    const updatePotParam = useCallback((key, value) => {
        setPotParams(prev => ({ ...prev, [key]: value }));
    }, []);

    const generatePlot = useCallback(async () => {
        setIsGenerating(true);
        setError(null);

        try {
            const response = await plotApi.generatePlotData(gridParams, plotType, potParams);
            setPlotData(response);
        } catch (err) {
            setError(err.message);
        } finally {
            setIsGenerating(false);
        }
    }, [gridParams, plotType, potParams]);

    const downloadPlot = useCallback(() => {
        const plotElement = document.querySelector('.js-plotly-plot');
        if (plotElement && window.Plotly) {
            window.Plotly.downloadImage(plotElement, {
                format: 'png',
                width: 1200,
                height: 800,
                filename: `plot_${Date.now()}`,
            });
        }
    }, []);

    // Memoized plot layout
    const plotLayout = useMemo(() => ({
        ...plotData?.layout,
        autosize: true,
        height: 600,
    }), [plotData]);

    const currentParamConfig = PARAM_CONFIGS[plotType];

    return (
        <div className="plot-generator">
            <div className="container">
                {/* Header */}
                <header className="header">
                    <h1 className="title">Interactive Plot Generator</h1>
                    <p className="subtitle">Generate plots for functions</p>
                </header>

                <div className="content-grid">
                    {/* Control Panel */}
                    <aside className="control-panel">
                        <h2 className="panel-title">
                            <BarChart3 className="icon" />
                            <span>Plot Settings</span>
                        </h2>

                        {/* Grid Parameters */}
                        <NumberInput
                            label="rMin"
                            value={gridParams.rMin}
                            onChange={(val) => updateGridParam('rMin', val)}
                        />
                        <NumberInput
                            label="rMax"
                            value={gridParams.rMax}
                            onChange={(val) => updateGridParam('rMax', val)}
                        />
                        <NumberInput
                            label="nGrid"
                            value={gridParams.nGrid}
                            onChange={(val) => updateGridParam('nGrid', val)}
                        />

                        {/* Plot Type Selector */}
                        <div className="input-group">
                            <label className="input-label">Plot Type</label>
                            <select
                                value={plotType}
                                onChange={(e) => setPlotType(e.target.value)}
                                className="select-field"
                            >
                                {PLOT_TYPES.map((type) => (
                                    <option key={type.value} value={type.value}>
                                        {type.label}
                                    </option>
                                ))}
                            </select>
                        </div>

                        {/* Dynamic Parameters */}
                        {currentParamConfig && (
                            <div className="params-section">
                                {currentParamConfig.map(({ key, label, step }) => (
                                    <NumberInput
                                        key={key}
                                        label={label}
                                        value={potParams[key]}
                                        step={step}
                                        onChange={(val) => updatePotParam(key, val)}
                                    />
                                ))}
                            </div>
                        )}

                        {/* Generate Button */}
                        <button
                            onClick={generatePlot}
                            disabled={isGenerating}
                            className="generate-btn"
                        >
                            {isGenerating ? (
                                <>
                                    <Loader className="icon animate-spin" />
                                    <span>Generating...</span>
                                </>
                            ) : (
                                <>
                                    <RefreshCw className="icon" />
                                    <span>Generate Plot</span>
                                </>
                            )}
                        </button>

                        {/* Error Display */}
                        {error && (
                            <div className="alert alert-error">
                                {error}
                            </div>
                        )}

                        {/* Info Box */}
                        <div className="info-box">
                            <p className="info-title">✨ Interactive Features:</p>
                            <ul className="info-list">
                                <li>• Zoom: Drag to select area</li>
                                <li>• Pan: Shift + drag</li>
                                <li>• Reset: Double click</li>
                                <li>• Hover for values</li>
                                <li>• Download via toolbar</li>
                            </ul>
                        </div>
                    </aside>

                    {/* Plot Display Area */}
                    <main className="plot-area">
                        <div className="plot-header">
                            <h2 className="plot-title">Visualization</h2>
                            {plotData && (
                                <button onClick={downloadPlot} className="download-btn">
                                    <Download className="icon" />
                                    <span>Download PNG</span>
                                </button>
                            )}
                        </div>

                        <div className="plot-container">
                            {isGenerating ? (
                                <div className="plot-placeholder">
                                    <div className="placeholder-content">
                                        <Loader className="placeholder-icon animate-spin" />
                                        <p className="placeholder-text">Generating plot...</p>
                                    </div>
                                </div>
                            ) : plotData ? (
                                <Plot
                                    data={plotData.data}
                                    layout={plotLayout}
                                    config={PLOT_CONFIG}
                                    style={{ width: '100%', height: '600px' }}
                                    useResizeHandler={true}
                                />
                            ) : (
                                <div className="plot-placeholder">
                                    <div className="placeholder-content">
                                        <BarChart3 className="placeholder-icon-large" />
                                        <p className="placeholder-text">
                                            Configure parameters and click "Generate Plot" to visualize
                                        </p>
                                    </div>
                                </div>
                            )}
                        </div>

                        {plotData && (
                            <div className="alert alert-success">
                                ✓ Plot generated successfully - Interact with the plot using your mouse!
                            </div>
                        )}
                    </main>
                </div>
            </div>
        </div>
    );
}