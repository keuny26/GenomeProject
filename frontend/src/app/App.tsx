import React, { useState, useEffect, useRef, useCallback } from 'react';
import ForceGraph2D, { ForceGraphMethods } from 'react-force-graph-2d';
import axios from 'axios';
import { Send, Upload, Database, X, Menu, Loader2, Trash2, RefreshCw, AlertCircle } from 'lucide-react';

// ⚠️ 백엔드 주소가 정확한지 다시 한번 확인하세요. (Docker 환경이라면 localhost 대신 IP가 필요할 수 있음)
const API_BASE_URL = "http://localhost:8000";

interface GraphNode {
  id: string;
  label: string;
  type: string;
  color: string;
  description: string;
  x?: number;
  y?: number;
}

interface GraphLink {
  source: string | any;
  target: string | any;
}

export default function App() {
  const [documents, setDocuments] = useState<{ id: string; name: string }[]>([]);
  const [graphData, setGraphData] = useState<{ nodes: GraphNode[], links: GraphLink[] }>({ nodes: [], links: [] });
  const [chatHistory, setChatHistory] = useState([{ role: 'ai', text: "시스템 준비 완료. PDF를 업로드해 주세요." }]);
  const [chatMessage, setChatMessage] = useState("");
  const [selectedNode, setSelectedNode] = useState<GraphNode | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [isSidebarOpen, setIsSidebarOpen] = useState(true);
  
  const chatEndRef = useRef<HTMLDivElement>(null);
  const fgRef = useRef<ForceGraphMethods | undefined>(undefined);

  const truncateText = (text: string, maxLength: number = 15) => {
    if (!text) return "";
    return text.length > maxLength ? text.substring(0, maxLength) + "..." : text;
  };

  const resetViewport = useCallback(() => {
    if (fgRef.current && graphData.nodes.length > 0) {
      fgRef.current.zoomToFit(800, 100);
      fgRef.current.centerAt(0, 0, 500);
    }
  }, [graphData.nodes]);

  useEffect(() => { chatEndRef.current?.scrollIntoView({ behavior: "smooth" }); }, [chatHistory]);

  useEffect(() => {
    if (fgRef.current && graphData.nodes.length > 0) {
      fgRef.current.d3Force('charge')?.strength(-500);
      fgRef.current.d3Force('link')?.distance((l: any) => 
        (l.source.type === 'doc' || l.target.type === 'doc') ? 180 : 60
      );
      fgRef.current.d3ReheatSimulation();
      // 데이터 로드 후 줌 자동 맞춤
      const timer = setTimeout(resetViewport, 500);
      return () => clearTimeout(timer);
    }
  }, [graphData, resetViewport]);

  // ✅ 파일 업로드 핸들러 보강
  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    // 초기화 및 로딩 시작
    setIsLoading(true);
    const formData = new FormData();
    formData.append('file', file);

    try {
      // 1. 서버 연결 테스트 (선택 사항이지만 진단에 도움됨)
      const res = await axios.post(`${API_BASE_URL}/upload`, formData, {
        headers: { 'Content-Type': 'multipart/form-data' },
        timeout: 30000 // 30초 타임아웃
      });

      console.log("Server Response:", res.data);

      if (res.data && res.data.graph) {
        const { nodes: newNodes, links: newLinks } = res.data.graph;

        // 그래프 데이터 업데이트
        setGraphData(prev => {
          const existingNodeIds = new Set(prev.nodes.map(n => n.id));
          const uniqueNewNodes = newNodes.filter((n: GraphNode) => !existingNodeIds.has(n.id));
          
          return {
            nodes: [...prev.nodes, ...uniqueNewNodes],
            links: [...prev.links, ...newLinks]
          };
        });

        setDocuments(prev => [...prev, { id: `doc_${Date.now()}`, name: file.name }]);
        setChatHistory(prev => [...prev, { role: 'ai', text: `✅ '${file.name}' 업로드 및 분석 성공!` }]);
      } else {
        throw new Error("서버에서 그래프 데이터를 받지 못했습니다.");
      }
    } catch (err: any) {
      console.error("Upload Error Details:", err);
      const errorMsg = err.response?.data?.detail || err.message || "서버 응답 없음";
      setChatHistory(prev => [...prev, { role: 'ai', text: `❌ 업로드 실패: ${errorMsg}` }]);
    } finally {
      setIsLoading(false);
      if (e.target) e.target.value = ''; // 입력창 초기화
    }
  };

  const handleSendMessage = async () => {
    if (!chatMessage.trim()) return;
    setChatHistory(prev => [...prev, { role: 'user', text: chatMessage }]);
    const currentMsg = chatMessage;
    setChatMessage("");
    try {
      const res = await axios.post(`${API_BASE_URL}/chat`, { query: currentMsg });
      setChatHistory(prev => [...prev, { role: 'ai', text: res.data.reply }]);
    } catch {
      setChatHistory(prev => [...prev, { role: 'ai', text: "질의에 실패했습니다." }]);
    }
  };

  return (
    <div className="flex h-screen w-full bg-white overflow-hidden font-sans text-slate-900">
      {isSidebarOpen && (
        <aside className="w-[350px] border-r flex flex-col bg-white shrink-0 shadow-lg z-20">
          <div className="p-6 border-b bg-slate-900 text-white flex items-center gap-3">
            <Database size={20} className="text-blue-400" />
            <h1 className="font-bold text-lg">GenomeGraph AI</h1>
          </div>

          <div className="flex-1 overflow-y-auto p-5 space-y-6">
            <label className={`flex flex-col items-center justify-center w-full h-32 border-2 border-dashed rounded-2xl cursor-pointer transition-all ${isLoading ? 'bg-slate-50 border-blue-400' : 'hover:bg-slate-50 border-slate-200'}`}>
              {isLoading ? (
                <div className="flex flex-col items-center">
                  <Loader2 className="animate-spin text-blue-500 mb-2" size={24} />
                  <span className="text-xs font-bold text-blue-500">서버 분석 중...</span>
                </div>
              ) : (
                <>
                  <Upload className="text-slate-400" size={24} />
                  <span className="text-xs font-bold mt-2 text-slate-500">PDF 업로드</span>
                </>
              )}
              <input type="file" className="hidden" onChange={handleFileUpload} accept=".pdf" disabled={isLoading} />
            </label>

            <section>
              <h3 className="text-[10px] font-black text-slate-400 uppercase tracking-widest mb-4">Analyzed Docs</h3>
              <div className="space-y-2">
                {documents.length === 0 && <p className="text-[10px] text-slate-400 text-center">분석된 문서가 없습니다.</p>}
                {documents.map(doc => (
                  <div key={doc.id} className="flex items-center justify-between p-3 bg-slate-50 border border-slate-100 rounded-xl">
                    <span className="text-xs font-semibold text-slate-600 truncate mr-2">📄 {doc.name}</span>
                    <button onClick={() => setDocuments(d => d.filter(x => x.id !== doc.id))} className="text-slate-300 hover:text-red-500">
                      <Trash2 size={14} />
                    </button>
                  </div>
                ))}
              </div>
            </section>
          </div>

          <div className="h-[380px] border-t bg-slate-50 flex flex-col">
            <div className="flex-1 overflow-y-auto p-4 space-y-3">
              {chatHistory.map((m, i) => (
                <div key={i} className={`flex ${m.role === 'user' ? 'justify-end' : 'justify-start'}`}>
                  <div className={`max-w-[85%] p-3 rounded-2xl text-[12px] shadow-sm ${m.role === 'user' ? 'bg-blue-600 text-white' : 'bg-white border text-slate-700'}`}>
                    {m.text}
                  </div>
                </div>
              ))}
              <div ref={chatEndRef} />
            </div>
            <div className="p-4 bg-white border-t flex gap-2">
              <input value={chatMessage} onChange={e => setChatMessage(e.target.value)} onKeyDown={e => e.key === 'Enter' && handleSendMessage()} className="flex-1 bg-slate-100 rounded-xl px-4 py-2 text-xs outline-none" placeholder="질문..." />
              <button onClick={handleSendMessage} className="p-2.5 bg-blue-600 text-white rounded-xl shadow-lg hover:bg-blue-700 transition-all"><Send size={16} /></button>
            </div>
          </div>
        </aside>
      )}

      <main className="flex-1 relative bg-[#F8FAFC] overflow-hidden">
        <div className="absolute top-6 left-6 z-30 flex gap-2">
          <button onClick={() => setIsSidebarOpen(!isSidebarOpen)} className="p-3 bg-white border rounded-2xl shadow-md hover:bg-slate-50"><Menu size={20} /></button>
          <button onClick={resetViewport} className="p-3 bg-white border rounded-2xl shadow-md hover:bg-slate-50 text-blue-600 flex items-center gap-2 font-bold text-xs"><RefreshCw size={16} /> 화면 맞춤</button>
        </div>

        {/* ✅ 데이터가 없을 때 표시할 안내 */}
        {graphData.nodes.length === 0 && !isLoading && (
          <div className="absolute inset-0 flex flex-col items-center justify-center text-slate-300 pointer-events-none">
            <Database size={64} className="mb-4 opacity-20" />
            <p className="text-sm font-medium">그래프 데이터가 없습니다.</p>
            <p className="text-xs">PDF를 업로드하면 지식 그래프가 생성됩니다.</p>
          </div>
        )}

        <ForceGraph2D
          ref={fgRef}
          graphData={graphData}
          nodeRelSize={1}
          nodeColor={n => (n as GraphNode).color}
          linkColor={l => (l.source.type === 'doc' || l.target.type === 'doc') ? "#E2E8F0" : "#94A3B8"}
          linkWidth={l => (l.source.type !== 'doc' && l.target.type !== 'doc') ? 3 : 1}
          linkDirectionalArrowLength={4}
          linkDirectionalArrowRelPos={1}
          onNodeClick={(node: any) => setSelectedNode(node)}
          nodeCanvasObject={(node, ctx, globalScale) => {
            const n = node as GraphNode;
            const isDoc = n.type === 'doc';
            const label = truncateText(n.label, isDoc ? 15 : 10);
            const fontSize = (isDoc ? 14 : 12) / globalScale;
            ctx.font = `${fontSize}px sans-serif`;
            ctx.fillStyle = n.color;
            ctx.beginPath(); 
            ctx.arc(node.x!, node.y!, isDoc ? 12 : 7, 0, 2 * Math.PI); 
            ctx.fill();
            
            if (selectedNode?.id === n.id) {
              ctx.strokeStyle = '#334155';
              ctx.lineWidth = 4 / globalScale;
              ctx.stroke();
            }
            ctx.fillStyle = isDoc ? "#0F172A" : "#475569"; 
            ctx.textAlign = 'center';
            ctx.fillText(label, node.x!, node.y! + (isDoc ? 24 : 18));
          }}
          cooldownTicks={100}
        />

        {selectedNode && (
          <div className="absolute bottom-10 right-10 w-80 bg-white p-6 rounded-[28px] shadow-2xl border z-50 transition-all">
            <div className="flex justify-between items-center mb-4">
              <span className={`px-2.5 py-1 rounded-lg text-[10px] font-black text-white ${selectedNode.type === 'gene' ? 'bg-blue-500' : selectedNode.type === 'disease' ? 'bg-red-500' : 'bg-slate-800'}`}>
                {selectedNode.type.toUpperCase()}
              </span>
              <button onClick={() => setSelectedNode(null)} className="p-1 hover:bg-slate-100 rounded-full">
                <X className="text-slate-400" size={18} />
              </button>
            </div>
            <h3 className="font-bold text-slate-900 text-xl mb-2">{selectedNode.label}</h3>
            <p className="text-[13px] text-slate-600 bg-slate-50 p-4 rounded-2xl leading-relaxed mb-6">{selectedNode.description}</p>
            <button onClick={() => window.open(`https://www.ncbi.nlm.nih.gov/search/all/?term=${selectedNode.label.replace(/🧬 |🏥 /g, '')}`)} className="w-full py-3.5 bg-slate-900 text-white rounded-2xl text-xs font-bold">상세 정보</button>
          </div>
        )}
      </main>
    </div>
  );
}