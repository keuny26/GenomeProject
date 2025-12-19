import { useState, useRef } from 'react';
import { motion, AnimatePresence } from 'motion/react';
import { 
  X, 
  Upload, 
  FileText, 
  Trash2, 
  Send,
  ChevronDown
} from 'lucide-react';

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
}

interface AnalysisItem {
  id: string;
  name: string;
  type: 'pdf' | 'document';
}

export function Sidebar({ isOpen, onClose }: SidebarProps) {
  const [gtexChecked, setGtexChecked] = useState(true);
  const [dbaseChecked, setDbaseChecked] = useState(true);
  const [filterOpen, setFilterOpen] = useState(true);
  const [chatMessage, setChatMessage] = useState('');
  const [analysisItems, setAnalysisItems] = useState<AnalysisItem[]>([
    { id: '1', name: 'Case Report Clinical impact of BRCA1 and BRCA2 in Breast...', type: 'pdf' }
  ]);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = e.target.files;
    if (files) {
      const newItems: AnalysisItem[] = Array.from(files).map((file, index) => ({
        id: `${Date.now()}-${index}`,
        name: file.name,
        type: 'pdf'
      }));
      setAnalysisItems([...analysisItems, ...newItems]);
    }
  };

  const handleDeleteItem = (id: string) => {
    setAnalysisItems(analysisItems.filter(item => item.id !== id));
  };

  const handleSendMessage = () => {
    if (chatMessage.trim()) {
      // 채팅 메시지 처리 로직
      console.log('메시지:', chatMessage);
      setChatMessage('');
    }
  };

  return (
    <AnimatePresence>
      {isOpen && (
        <motion.aside
          initial={{ x: -320 }}
          animate={{ x: 0 }}
          exit={{ x: -320 }}
          transition={{ type: 'spring', damping: 25, stiffness: 200 }}
          className="w-80 bg-white border-r border-gray-200 flex flex-col h-full shadow-lg z-10 fixed md:relative left-0 top-0"
        >
          {/* 헤더 */}
          <div className="bg-blue-600 text-white p-4 flex items-center justify-between">
            <motion.div 
              initial={{ opacity: 0, x: -20 }}
              animate={{ opacity: 1, x: 0 }}
              transition={{ delay: 0.2 }}
              className="flex items-center gap-2"
            >
              <div className="w-6 h-6 bg-white rounded flex items-center justify-center">
                <span className="text-blue-600 font-bold text-sm">G</span>
              </div>
              <h1 className="font-semibold">Genome AI</h1>
            </motion.div>
            <button
              onClick={onClose}
              className="p-1 hover:bg-blue-700 rounded transition-colors"
              aria-label="사이드바 닫기"
            >
              <X className="w-5 h-5" />
            </button>
          </div>

          {/* PDF 업로드 버튼 */}
          <div className="p-4 border-b border-gray-200">
            <input
              ref={fileInputRef}
              type="file"
              accept=".pdf,application/pdf"
              multiple
              onChange={handleFileUpload}
              className="hidden"
            />
            <button
              onClick={() => fileInputRef.current?.click()}
              className="w-full bg-blue-600 hover:bg-blue-700 text-white py-2.5 px-4 rounded-lg flex items-center justify-center gap-2 transition-colors"
            >
              <Upload className="w-4 h-4" />
              <span>PDF 분석 추가</span>
            </button>
          </div>

          {/* 데이터 필터 */}
          <div className="border-b border-gray-200">
            <button
              onClick={() => setFilterOpen(!filterOpen)}
              className="w-full px-4 py-3 flex items-center justify-between hover:bg-gray-50 transition-colors"
            >
              <span className="text-sm font-medium">데이터 필터</span>
              <ChevronDown 
                className={`w-4 h-4 transition-transform ${filterOpen ? 'rotate-180' : ''}`}
              />
            </button>
            
            <AnimatePresence>
              {filterOpen && (
                <motion.div
                  initial={{ height: 0, opacity: 0 }}
                  animate={{ height: 'auto', opacity: 1 }}
                  exit={{ height: 0, opacity: 0 }}
                  transition={{ duration: 0.2 }}
                  className="overflow-hidden"
                >
                  <div className="px-4 pb-3 space-y-2">
                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={gtexChecked}
                        onChange={(e) => setGtexChecked(e.target.checked)}
                        className="w-4 h-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm">GTRF</span>
                    </label>
                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={dbaseChecked}
                        onChange={(e) => setDbaseChecked(e.target.checked)}
                        className="w-4 h-4 rounded border-gray-300 text-blue-600 focus:ring-blue-500"
                      />
                      <span className="text-sm">DBASE</span>
                    </label>
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>

          {/* 분석 항목 */}
          <div className="flex-1 overflow-y-auto">
            <div className="px-4 py-3">
              <h3 className="text-sm font-medium text-gray-700 mb-2">분석 항목</h3>
              <div className="space-y-2">
                {analysisItems.map((item) => (
                  <div
                    key={item.id}
                    className="group flex items-start gap-2 p-2 rounded-lg hover:bg-gray-50 transition-colors"
                  >
                    <FileText className="w-4 h-4 text-blue-600 mt-0.5 flex-shrink-0" />
                    <span className="flex-1 text-xs text-gray-700 line-clamp-2">
                      {item.name}
                    </span>
                    <button
                      onClick={() => handleDeleteItem(item.id)}
                      className="opacity-0 group-hover:opacity-100 p-1 hover:bg-red-50 rounded transition-all"
                      aria-label="삭제"
                    >
                      <Trash2 className="w-3.5 h-3.5 text-red-500" />
                    </button>
                  </div>
                ))}
              </div>
              {analysisItems.length === 0 && (
                <p className="text-xs text-gray-400 text-center py-4">
                  아직 분석할 항목이 없습니다
                </p>
              )}
            </div>
          </div>

          {/* 채팅 입력 */}
          <div className="border-t border-gray-200 p-4">
            <div className="text-xs text-gray-500 mb-2">
              궁금하신 내용을 자유롭게 입력해주세요
            </div>
            <div className="flex gap-2">
              <input
                type="text"
                value={chatMessage}
                onChange={(e) => setChatMessage(e.target.value)}
                onKeyDown={(e) => {
                  if (e.key === 'Enter' && !e.shiftKey) {
                    e.preventDefault();
                    handleSendMessage();
                  }
                }}
                placeholder="질문하기..."
                className="flex-1 px-3 py-2 border border-gray-300 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
              />
              <button
                onClick={handleSendMessage}
                className="p-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 transition-colors disabled:opacity-50 disabled:cursor-not-allowed"
                disabled={!chatMessage.trim()}
                aria-label="전송"
              >
                <Send className="w-4 h-4" />
              </button>
            </div>
          </div>
        </motion.aside>
      )}
    </AnimatePresence>
  );
}