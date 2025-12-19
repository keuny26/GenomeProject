import { useEffect, useRef, useState } from 'react';
import { motion } from 'motion/react';

interface Node {
  id: string;
  label: string;
  x: number;
  y: number;
  type: 'center' | 'gene' | 'document';
  vx?: number;
  vy?: number;
  isDragging?: boolean;
}

interface Edge {
  source: string;
  target: string;
}

export function NetworkGraph() {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [nodes, setNodes] = useState<Node[]>([]);
  const [edges, setEdges] = useState<Edge[]>([]);
  const [hoveredNode, setHoveredNode] = useState<string | null>(null);
  const [draggedNode, setDraggedNode] = useState<string | null>(null);
  const [isDraggingCanvas, setIsDraggingCanvas] = useState(false);
  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const animationRef = useRef<number>();
  const mousePos = useRef({ x: 0, y: 0 });
  const dragStart = useRef({ x: 0, y: 0 });
  const lastOffset = useRef({ x: 0, y: 0 });

  // 초기 네트워크 데이터 생성
  useEffect(() => {
    const centerNode: Node = {
      id: 'center',
      label: 'TP53',
      x: 400,
      y: 300,
      type: 'center'
    };

    const geneNames = ['BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'RAD51', 'MDM2', 'CDKN1A', 'BAX', 'PUMA', 'NOXA'];
    const geneNodes: Node[] = geneNames.map((name, i) => {
      const angle = (i / geneNames.length) * Math.PI * 2;
      const radius = 150;
      return {
        id: `gene-${i}`,
        label: name,
        x: 400 + Math.cos(angle) * radius,
        y: 300 + Math.sin(angle) * radius,
        type: 'gene'
      };
    });

    const documentNodes: Node[] = [
      { id: 'doc1', label: 'COVID-19 & CANC...', x: 700, y: 200, type: 'document' },
      { id: 'doc2', label: 'ADMIN FINANC...', x: 650, y: 250, type: 'document' },
      { id: 'doc3', label: 'Case Report ...', x: 720, y: 300, type: 'document' },
      { id: 'doc4', label: 'CANCER', x: 780, y: 240, type: 'document' }
    ];

    const centerNode2: Node = {
      id: 'center2',
      label: 'BRCA',
      x: 850,
      y: 320,
      type: 'center'
    };

    const cluster2Genes = ['PIK3CA', 'PTEN', 'AKT1', 'MTOR', 'TSC1', 'TSC2'];
    const cluster2Nodes: Node[] = cluster2Genes.map((name, i) => {
      const angle = (i / cluster2Genes.length) * Math.PI * 2;
      const radius = 100;
      return {
        id: `gene2-${i}`,
        label: name,
        x: 850 + Math.cos(angle) * radius,
        y: 320 + Math.sin(angle) * radius,
        type: 'gene'
      };
    });

    const allNodes = [centerNode, ...geneNodes, ...documentNodes, centerNode2, ...cluster2Nodes];
    
    const geneEdges: Edge[] = geneNodes.map(node => ({
      source: 'center',
      target: node.id
    }));

    const docEdges: Edge[] = documentNodes.map(node => ({
      source: 'center2',
      target: node.id
    }));

    const cluster2Edges: Edge[] = cluster2Nodes.map(node => ({
      source: 'center2',
      target: node.id
    }));

    setNodes(allNodes);
    setEdges([...geneEdges, ...docEdges, ...cluster2Edges]);
  }, []);

  // 캔버스 그리기
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const resizeCanvas = () => {
      const rect = canvas.getBoundingClientRect();
      canvas.width = rect.width * window.devicePixelRatio;
      canvas.height = rect.height * window.devicePixelRatio;
      ctx.scale(window.devicePixelRatio, window.devicePixelRatio);
    };

    resizeCanvas();
    window.addEventListener('resize', resizeCanvas);

    const draw = () => {
      const rect = canvas.getBoundingClientRect();
      ctx.clearRect(0, 0, rect.width, rect.height);

      // offset 적용
      ctx.save();
      ctx.translate(offset.x, offset.y);

      // 엣지 그리기
      edges.forEach(edge => {
        const sourceNode = nodes.find(n => n.id === edge.source);
        const targetNode = nodes.find(n => n.id === edge.target);
        
        if (sourceNode && targetNode) {
          ctx.beginPath();
          ctx.moveTo(sourceNode.x, sourceNode.y);
          ctx.lineTo(targetNode.x, targetNode.y);
          ctx.strokeStyle = '#FFA500';
          ctx.lineWidth = 2;
          ctx.stroke();
        }
      });

      // 노드 그리기
      nodes.forEach(node => {
        const isHovered = hoveredNode === node.id;
        
        // 노드 원
        ctx.beginPath();
        ctx.arc(node.x, node.y, isHovered ? 10 : 8, 0, Math.PI * 2);
        
        if (node.type === 'center') {
          ctx.fillStyle = '#FF6B6B';
        } else if (node.type === 'document') {
          ctx.fillStyle = '#FF4444';
        } else {
          ctx.fillStyle = '#4A90E2';
        }
        ctx.fill();
        
        // 외곽선
        ctx.strokeStyle = 'white';
        ctx.lineWidth = 2;
        ctx.stroke();

        // 라벨
        ctx.fillStyle = '#333';
        ctx.font = '12px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'top';
        ctx.fillText(node.label, node.x, node.y + 12);
      });

      ctx.restore();

      animationRef.current = requestAnimationFrame(draw);
    };

    draw();

    return () => {
      window.removeEventListener('resize', resizeCanvas);
      if (animationRef.current) {
        cancelAnimationFrame(animationRef.current);
      }
    };
  }, [nodes, edges, hoveredNode, offset]);

  // 마우스 인터랙션
  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    // 캔버스 드래그 중일 때
    if (isDraggingCanvas) {
      const dx = mouseX - dragStart.current.x;
      const dy = mouseY - dragStart.current.y;
      setOffset({ x: lastOffset.current.x + dx, y: lastOffset.current.y + dy });
      return;
    }

    // 노드 드래그 중일 때 (offset 고려)
    if (draggedNode) {
      const worldX = mouseX - offset.x;
      const worldY = mouseY - offset.y;
      setNodes(prevNodes => 
        prevNodes.map(node => 
          node.id === draggedNode 
            ? { ...node, x: worldX, y: worldY }
            : node
        )
      );
      mousePos.current = { x: worldX, y: worldY };
      return;
    }

    // 호버 감지 (offset 고려)
    const worldX = mouseX - offset.x;
    const worldY = mouseY - offset.y;
    
    let found = false;
    for (const node of nodes) {
      const dx = worldX - node.x;
      const dy = worldY - node.y;
      const distance = Math.sqrt(dx * dx + dy * dy);
      
      if (distance < 15) {
        setHoveredNode(node.id);
        found = true;
        break;
      }
    }

    if (!found) {
      setHoveredNode(null);
    }

    mousePos.current = { x: worldX, y: worldY };
  };

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const canvas = canvasRef.current;
    if (!canvas) return;

    const rect = canvas.getBoundingClientRect();
    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    // offset 고려한 월드 좌표
    const worldX = mouseX - offset.x;
    const worldY = mouseY - offset.y;

    // 노드 클릭 감지
    let nodeClicked = false;
    for (const node of nodes) {
      const dx = worldX - node.x;
      const dy = worldY - node.y;
      const distance = Math.sqrt(dx * dx + dy * dy);
      
      if (distance < 15) {
        setDraggedNode(node.id);
        nodeClicked = true;
        break;
      }
    }

    // 노드가 아닌 빈 공간 클릭 시 캔버스 드래그 시작
    if (!nodeClicked) {
      setIsDraggingCanvas(true);
      dragStart.current = { x: mouseX, y: mouseY };
      lastOffset.current = offset;
    }
  };

  const handleMouseUp = () => {
    setDraggedNode(null);
    setIsDraggingCanvas(false);
    lastOffset.current = offset;
  };

  const handleMouseLeave = () => {
    setHoveredNode(null);
    setDraggedNode(null);
    setIsDraggingCanvas(false);
  };

  return (
    <div className="size-full relative bg-gradient-to-br from-gray-50 to-blue-50">
      <canvas
        ref={canvasRef}
        onMouseMove={handleMouseMove}
        onMouseDown={handleMouseDown}
        onMouseUp={handleMouseUp}
        onMouseLeave={handleMouseLeave}
        className="size-full cursor-pointer"
      />
      
      {/* 범례 */}
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ delay: 0.5 }}
        className="absolute bottom-6 left-6 bg-white rounded-lg shadow-md p-4 space-y-2"
      >
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-[#FF6B6B]"></div>
          <span className="text-xs text-gray-700">중심 유전자</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-[#4A90E2]"></div>
          <span className="text-xs text-gray-700">관련 유전자</span>
        </div>
        <div className="flex items-center gap-2">
          <div className="w-3 h-3 rounded-full bg-[#FF4444]"></div>
          <span className="text-xs text-gray-700">문서</span>
        </div>
      </motion.div>
    </div>
  );
}