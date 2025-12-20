import React, { useState, useMemo, useCallback } from 'react';
import * as math from 'mathjs';

const presetFunctions = [
  { label: 'Quadratic', expr: 'x^2' },
  { label: 'Sine Wave', expr: 'sin(x)' },
  { label: 'Cosine Wave', expr: 'cos(x)' },
  { label: 'Exponential', expr: 'e^x' },
  { label: 'Logarithm', expr: 'log(x)' },
  { label: 'Cubic', expr: 'x^3 - 3x' },
  { label: 'Tangent', expr: 'tan(x)' },
  { label: 'Square Root', expr: 'sqrt(x)' },
  { label: 'Gaussian', expr: 'e^(-x^2)' },
  { label: 'Sinc', expr: 'sin(x)/x' },
];

const mathHelp = [
  { category: 'Basic', ops: ['+', '-', '*', '/', '^'] },
  { category: 'Trig', ops: ['sin', 'cos', 'tan', 'asin', 'acos', 'atan'] },
  { category: 'Other', ops: ['sqrt', 'log', 'exp', 'abs', 'floor', 'ceil'] },
  { category: 'Constants', ops: ['pi', 'e'] },
];

export default function FunctionVisualizer() {
  const [expression, setExpression] = useState('sin(x) * x');
  const [xMin, setXMin] = useState(-10);
  const [xMax, setXMax] = useState(10);
  const [nPoints, setNpoints] = useState(200);
  const [error, setError] = useState(null);

  const parseExpression = useCallback((expr) => {
    try {
      // Pre-process the expression for common mathematical conventions
      let processed = expr
        .replace(/(\d)([a-zA-Z])/g, '$1*$2')  // 2x -> 2*x
        .replace(/([a-zA-Z])(\d)/g, '$1*$2')  // x2 -> x*2 (less common but handle it)
        .replace(/\)\(/g, ')*(')              // )(  -> )*(
        .replace(/(\d)\(/g, '$1*(')           // 2(  -> 2*(
        .replace(/\)(\d)/g, ')*$1')           // )2  -> )*2
        .replace(/([a-zA-Z])\(/g, (match, p1) => {
          // Don't add * for known functions
          const funcs = ['sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'sinh', 'cosh', 'tanh', 
                        'sqrt', 'log', 'ln', 'exp', 'abs', 'floor', 'ceil', 'round', 'sign'];
          if (funcs.some(f => expr.includes(f + '('))) {
            return match;
          }
          return p1 + '*(';
        });
      
      const compiled = math.compile(processed);
      return { compiled, processed, error: null };
    } catch (e) {
      return { compiled: null, processed: expr, error: e.message };
    }
  }, []);

  const { data, parsedExpr, parseError } = useMemo(() => {
    const { compiled, processed, error } = parseExpression(expression);
    
    if (error) {
      setError(error);
      return { data: [], parsedExpr: processed, parseError: error };
    }
    
    setError(null);
    const points = [];
    const step = (xMax - xMin) / nPoints;
    
    for (let i = 0; i <= nPoints; i++) {
      const x = xMin + i * step;
      try {
        let y = compiled.evaluate({ x });
        
        // Handle special cases
        if (!isFinite(y) || isNaN(y)) {
          y = null;
        } else if (Math.abs(y) > 1e10) {
          y = null; // Clip extremely large values
        }
        
        points.push({ x: parseFloat(x.toFixed(4)), y });
      } catch {
        points.push({ x: parseFloat(x.toFixed(4)), y: null });
      }
    }
    
    return { data: points, parsedExpr: processed, parseError: null };
  }, [expression, xMin, xMax, nPoints, parseExpression]);

  const yDomain = useMemo(() => {
    const validYs = data.filter(p => p.y !== null).map(p => p.y);
    if (validYs.length === 0) return [-10, 10];
    
    const minY = Math.min(...validYs);
    const maxY = Math.max(...validYs);
    const padding = (maxY - minY) * 0.1 || 1;
    
    return [
      Math.max(minY - padding, -100),
      Math.min(maxY + padding, 100)
    ];
  }, [data]);

  return (
    <div style={{
        padding: '24px',
    }}>
        <div style={{
            backgroundSize: '40px 40px',
            pointerEvents: 'auto',
            position: 'fixed',
            inset: 0,

        }}>
            <h1 style={{ textAlign: 'center', marginBottom: '32px', padding: '20px'}}>Custom Function Visualizer</h1>
            <label style={{ 
                paddingLeft: '30px',
                display: 'block', 
                fontSize: '1.5rem', 
                color: '#63b3ed',
                letterSpacing: '0.15em',
                textTransform: 'uppercase',
                marginBottom: '8px',
            }}>
                Function Expression
            </label>

            <div style={{ display: 'flex', gap: '12px', alignItems: 'center', padding: '30px' }}>
                <span style={{ 
                    fontSize: '1.5rem', 
                    color: '#dbe72ee0',
                    fontStyle: 'italic',
                }}>
                    f(x) =
                </span>
            
                <input
                    type="text"
                    value={expression}
                    onChange={(e) => setExpression(e.target.value)}
                    placeholder="e.g., sin(x) * x^2"
                    style={{
                        flex: 1,
                        background: 'rgba(0, 0, 0, 0.4)',
                        border: error ? '2px solid #fc8181' : '2px solid rgba(99, 179, 237, 0.3)',
                        borderRadius: '8px',
                        padding: '16px 20px',
                        fontSize: '1.25rem',
                        color: '#fff',
                        fontFamily: 'inherit',
                        outline: 'none',
                        transition: 'border-color 0.2s, box-shadow 0.2s',
                    }}
                    onFocus={(e) => {
                        e.target.style.borderColor = '#63b3ed';
                        e.target.style.boxShadow = '0 0 20px rgba(99, 179, 237, 0.3)';
                    }}
                    onBlur={(e) => {
                        e.target.style.borderColor = error ? '#fc8181' : 'rgba(99, 179, 237, 0.3)';
                        e.target.style.boxShadow = 'none';
                    }}
                />
            </div>
            {error && (
                <div style={{
                    marginTop: '12px',
                    padding: '12px 16px',
                    background: 'rgba(252, 129, 129, 0.1)',
                    border: '1px solid rgba(252, 129, 129, 0.3)',
                    borderRadius: '8px',
                    color: '#fc8181',
                    fontSize: '0.875rem',
                }}>
                    ⚠ Parse Error: {error}
                </div>
            )}

            {/* Preset Functions */}
            <div style={{ marginTop: '0px', paddingLeft: '30px' }}>
                <span style={{ 
                    fontSize: '0.7rem', 
                    color: '#718096',
                    letterSpacing: '0.1em',
                    textTransform: 'uppercase',
                }}>
                    Quick Presets:
                </span>
                <div style={{ 
                    display: 'flex',  
                    gap: '20px', 
                    marginTop: '8px' 
                }}>
                    {presetFunctions.map((preset) => (
                        <button
                            key={preset.expr}
                            onClick={() => setExpression(preset.expr)}
                            style={{
                                border: '1px solid rgba(99, 179, 237, 0.3)',
                                borderRadius: '20px',
                                padding: '6px 14px',
                                color: expression === preset.expr ? '#fff' : '#a0aec0',
                                fontSize: '0.75rem',
                                cursor: 'pointer',
                                transition: 'all 0.2s',
                                fontFamily: 'inherit',
                            }}
                            onMouseEnter={(e) => {
                                if (expression !== preset.expr) {
                                    e.target.style.background = 'rgba(99, 179, 237, 0.2)';
                                    e.target.style.color = '#fff';
                                }
                            }}
                            onMouseLeave={(e) => {
                                if (expression !== preset.expr) {
                                    e.target.style.background = 'rgba(255, 255, 255, 0.05)';
                                    e.target.style.color = '#a0aec0';
                                }
                            }}
                        >
                            {preset.label}
                        </button>
                    ))}
                </div>
            </div>

            {/* Controls Row */}
            <div style={{
                display: 'flex',
                flexDirection: 'column',
                gap: '16px',
                marginBottom: '24px',
                alignItems: 'flex-start',
                padding: '30px'
            }}>
                {/* X Range Controls */}
                <div style={{
                    background: 'rgba(26, 32, 44, 0.6)',
                    borderRadius: '12px',
                    border: '1px solid rgba(99, 179, 237, 0.15)',
                    padding: '16px',
                    alignContent: 'center',
                }}>
                    <label style={{ 
                        fontSize: '0.7rem', 
                        color: '#63b3ed',
                        letterSpacing: '0.1em',
                        textTransform: 'uppercase',
                    }}>
                        X Range
                    </label>
                    <div style={{ display: 'flex', gap: '12px', marginTop: '8px', alignItems: 'center', paddingBottom: '8px'}}>
                        <input
                            type="number"
                            value={xMin}
                            onChange={(e) => setXMin(parseFloat(e.target.value) || -10)}
                            style={{
                                width: '80px',
                                background: 'rgba(0, 0, 0, 0.4)',
                                border: '1px solid rgba(99, 179, 237, 0.3)',
                                borderRadius: '6px',
                                padding: '8px 12px',
                                color: '#fff',
                                fontSize: '0.875rem',
                                fontFamily: 'inherit',
                                textAlign: 'center',
                            }}
                        />
                        <span style={{ color: '#fff' }}>to</span>
                        <input
                            type="number"
                            value={xMax}
                            onChange={(e) => setXMax(parseFloat(e.target.value) || 10)}
                            style={{
                                width: '80px',
                                background: 'rgba(0, 0, 0, 0.4)',
                                border: '1px solid rgba(99, 179, 237, 0.3)',
                                borderRadius: '6px',
                                padding: '8px 12px',
                                color: '#fff',
                                fontSize: '0.875rem',
                                fontFamily: 'inherit',
                                textAlign: 'center',
                            }}
                        />
                    </div>

                    <label style={{ 
                        paddingTop: '20px',
                        fontSize: '0.7rem', 
                        color: '#63b3ed',
                        letterSpacing: '0.1em',
                        textTransform: 'uppercase',
                    }}>
                        Resolution: {nPoints} points
                    </label>
                    <div style={{ display: 'flex', gap: '12px', marginTop: '8px', alignItems: 'center' }}>
                        <input
                            type="number"
                            value={nPoints}
                            onChange={(e) => setNPoints(parseInt(e.target.value) || 200)}
                            style={{
                                width: '80px',
                                background: 'rgba(0, 0, 0, 0.4)',
                                border: '1px solid rgba(99, 179, 237, 0.3)',
                                borderRadius: '6px',
                                padding: '8px 12px',
                                color: '#fff',
                                fontSize: '0.875rem',
                                fontFamily: 'inherit',
                                textAlign: 'center',
                            }}
                        />
                    </div>
                </div>

                {/* Math Reference */}
                <div style={{
                    background: 'rgba(26, 32, 44, 0.6)',
                    borderRadius: '12px',
                    border: '1px solid rgba(99, 179, 237, 0.15)',
                    padding: '16px',
                    gridColumn: 'span 1',
                }}>
                    <label style={{ 
                        fontSize: '1rem', 
                        color: '#63b3ed',
                        letterSpacing: '0.2em',
                        textTransform: 'uppercase',
                    }}>
                        Operators & Functions
                    </label>
                    <div style={{ marginTop: '8px', fontSize: '1rem' }}>
                        {mathHelp.map((cat) => (
                            <div key={cat.category} style={{ marginBottom: '4px' }}>
                                <span style={{ color: '#718096' }}>{cat.category}: </span>
                                <span style={{ color: '#e2e8f0' }}>{cat.ops.join(', ')}</span>
                            </div>
                        ))}
                    </div>
                </div>
            </div>
        </div>    
    </div>
        
  );
}