import { BrowserRouter, Routes, Route } from 'react-router-dom'
import Homepage from './pages/HomePage'
import PlotGenerator from './components/PlotGenerator.jsx'
import CalculatorView from "./calculator/calculatorView.jsx";
import './index.css'

function App() {
    return (
        <BrowserRouter>
            <Routes>
                <Route path="/" element={<Homepage />} />
                <Route path="/plots" element={<PlotGenerator />} />
                <Route path="/Calculator" element={<CalculatorView />} />
                <Route path="/webmanual" />
            </Routes>
        </BrowserRouter>
    )
}

export default App
